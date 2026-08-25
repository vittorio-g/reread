/* app.js — UI logic for the reread web tool. All local, no network.
 * (JS namespace stays `ReReRe` for back-compat; user-facing name is reread.) */
(function () {
  "use strict";
  const $ = id => document.getElementById(id);
  const drop = $("drop"), fileInput = $("file");
  let lastResult = null;      // {res, ids, header, excluded, fileName}
  let worker = null;

  /* ---- bespoke line-art icons (no emoji anywhere) ---- */
  const ICO = {
    check: '<svg class="ico" viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"><path d="M20 6.5 9.5 17 4 11.5"/></svg>',
    warn:  '<svg class="ico" viewBox="0 0 24 24" width="15" height="15" fill="none" stroke="currentColor" stroke-width="2.1" stroke-linecap="round" stroke-linejoin="round"><path d="M12 3.5 21.5 20H2.5z"/><path d="M12 10v4.5"/><path d="M12 17.6h.01"/></svg>',
    flag:  '<svg class="ico" viewBox="0 0 24 24" width="12" height="12" fill="none" stroke="currentColor" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"><path d="M5.5 21V3.5"/><path d="M5.5 4.2h11l-2.2 3.4 2.2 3.4h-11"/></svg>'
  };

  /* ---------------- helpers ---------------- */
  function setError(msg) { $("error").textContent = msg || ""; }
  function setProgress(p, label) {
    $("progressWrap").style.display = p == null ? "none" : "block";
    if (p != null) {
      $("progressFill").style.width = Math.round(p * 100) + "%";
      $("progressLabel").textContent = label || "Computing…";
    }
  }
  /* manual expected-rate override (Advanced). Blank/invalid -> null = automatic. */
  function manualRate() {
    const pct = parseFloat($("pPrev").value);
    return (Number.isFinite(pct) && pct > 0 && pct < 100) ? pct / 100 : null;
  }
  /* sensitivity: one shipped cut ("standard", z=2.5); a raw numeric z is the
   * power-user override for the out-of-envelope case. */
  function sensSetting() {
    const cz = $("pSensZ") ? parseFloat($("pSensZ").value) : NaN;
    return Number.isFinite(cz) ? cz : "standard";
  }
  /* delimiter: a free-typed single character wins over the preset dropdown. */
  function delimSetting() {
    const c = $("pDelimCustom") ? $("pDelimCustom").value : "";
    if (c && c.length) return c.charAt(0);
    return $("pDelim").value || null;
  }
  function opts() {
    const o = {
      corProp:    parseFloat($("pCorProp").value) || 0.03,
      iterations: parseInt($("pIter").value, 10) || 200,
      minPairs:   parseInt($("pMinPairs").value, 10) || 15,
      seed:       parseInt($("pSeed").value, 10) || 1,
      delimiter:  delimSetting(),
      sensitivity: sensSetting()
    };
    const mr = manualRate();
    if (mr != null) o.prevalence = mr;      // else: automatic attentive-mode calibration
    return o;
  }
  /* re-derive flag cutoff & flags from the ensemble probability + a target
   * careless rate (top-share). No recompute needed when the rate changes. */
  function applyPrevalence(res, prev) {
    const p = res.p, N = p.length;
    const srt = [...p].sort((a, b) => a - b);
    const cut = Math.min(N - 1, Math.max(0, Math.floor((1 - prev) * N)));
    const thr = srt[cut];
    res.threshold = thr;
    res.flagged = p.map(v => v >= thr);
    res.nFlagged = res.flagged.filter(Boolean).length;
    res.params.prevalence = prev;
    return thr;
  }

  /* ---------------- file intake ---------------- */
  drop.addEventListener("click", () => fileInput.click());
  drop.addEventListener("dragover", e => { e.preventDefault(); drop.classList.add("hover"); });
  drop.addEventListener("dragleave", () => drop.classList.remove("hover"));
  drop.addEventListener("drop", e => {
    e.preventDefault(); drop.classList.remove("hover");
    if (e.dataTransfer.files.length) handleFile(e.dataTransfer.files[0]);
  });
  fileInput.addEventListener("change", () => {
    if (fileInput.files.length) handleFile(fileInput.files[0]);
  });

  function handleFile(f) {
    setError("");
    if (/\.(sav|ods)$/i.test(f.name)) {
      setError("SPSS/ODS files are not supported — please export as CSV or XLSX and drop that.");
      return;
    }
    if (/\.(xlsx|xls|xlsm|xlsb)$/i.test(f.name)) { handleExcel(f); return; }
    const reader = new FileReader();
    reader.onload = () => {
      $("fileInfo").textContent = f.name + " (" + Math.round(f.size / 1024) + " KB)";
      compute(reader.result, f.name);
    };
    reader.onerror = () => setError("Could not read the file.");
    reader.readAsText(f);
  }

  /* Excel: read locally with SheetJS, ask which sheet, convert to CSV. */
  function handleExcel(f) {
    if (typeof XLSX === "undefined") {
      setError("The Excel reader failed to load. Please refresh, or export the file as CSV.");
      return;
    }
    const reader = new FileReader();
    reader.onerror = () => setError("Could not read the file.");
    reader.onload = () => {
      let wb;
      try { wb = XLSX.read(reader.result, { type: "array" }); }
      catch (e) { setError("Could not parse this Excel file (" + (e.message || e) + "). Try re-saving it, or export as CSV."); return; }
      const sheets = (wb.SheetNames || []).filter(name =>
        wb.Sheets[name] && wb.Sheets[name]["!ref"]);   // skip truly empty sheets
      if (!sheets.length) { setError("This workbook has no non-empty sheets."); return; }
      const pick = name => {
        const csv = XLSX.utils.sheet_to_csv(wb.Sheets[name], { FS: "," });
        $("fileInfo").textContent = f.name + " › sheet “" + name + "” (" + Math.round(f.size / 1024) + " KB)";
        compute(csv, f.name.replace(/\.\w+$/, "") + "_" + name + ".csv");
      };
      if (sheets.length === 1) pick(sheets[0]);
      else showSheetPicker(sheets, wb, pick);
    };
    reader.readAsArrayBuffer(f);
  }

  function showSheetPicker(sheets, wb, onPick) {
    const modal = $("sheetModal"), list = $("sheetList");
    list.innerHTML = "";
    sheets.forEach(name => {
      const ws = wb.Sheets[name];
      let rows = 0, cols = 0;
      try {
        const rg = XLSX.utils.decode_range(ws["!ref"]);
        rows = rg.e.r - rg.s.r + 1; cols = rg.e.c - rg.s.c + 1;
      } catch (e) {}
      const btn = document.createElement("button");
      btn.className = "ghost";
      btn.style.cssText = "text-align:left;padding:12px 14px;width:100%";
      btn.innerHTML = "<b>" + escapeHtml(name) + "</b> <span style='color:var(--slate);font-weight:400'>— " +
        rows + " rows × " + cols + " cols</span>";
      btn.addEventListener("click", () => { hideSheetPicker(); onPick(name); });
      list.appendChild(btn);
    });
    modal.style.display = "flex";
  }
  function hideSheetPicker() { $("sheetModal").style.display = "none"; }
  $("sheetCancel").addEventListener("click", hideSheetPicker);
  $("sheetModal").addEventListener("click", e => { if (e.target === $("sheetModal")) hideSheetPicker(); });

  /* ---------------- demo ----------------
   * The reread validation sample (Study 1): a real five-instrument battery
   * answered by real people, some of whom were instructed to answer carelessly
   * on a known share of it. Nothing is injected or simulated — the ground truth
   * was collected, and it travels in the row id.
   *
   * demo_data.js carries the engine's OUTPUT only, computed offline at these
   * same default settings (see webapp/extract_demo_study1.js); the responses are
   * not published. So the demo skips the compute step and renders the stored
   * result — everything downstream (plot, table, downloads, the live re-flag
   * controls, which work off eta) behaves exactly as it does on your own file. */
  $("demoBtn").addEventListener("click", () => {
    setError("");
    const D = self.RERERE_DEMO_RESULT, m = self.RERERE_DEMO_META || {};
    if (!D || !D.res) { setError("Demo data failed to load — please refresh."); return; }
    const scales = m.scales ? Object.entries(m.scales).map(([k, v]) => v + " " + k).join(", ") : "";
    const rate = Math.round(100 * m.careless / m.n);
    $("fileInfo").textContent = "demo: the reread validation sample — " + m.n +
      " real respondents × " + m.items + " items (" + scales + "). " + m.careful +
      " answered honestly; " + m.careless + " were instructed to answer carelessly on part of the " +
      "battery (" + m.byGroup["careless_50"] + " on 50%, " + m.byGroup["careless_75"] + " on 75%, " +
      m.byGroup["careless_100"] + " on 100%). The row id is the ground truth — P = careful, " +
      "C50/C75/C100 = instructed careless — so you can check every flag against it. Careless " +
      "respondents were deliberately over-sampled here (" + rate + "% of the sample, against " +
      "8–12% in a typical survey), because the study needed enough of them to fit the method. " +
      "These results were computed on the real responses at the default settings; the responses " +
      "themselves are not published here, so this view is the tool's output, not its input.";
    setProgress(null);
    renderAll(D, "reread_demo.csv");
  });

  /* ---------------- compute via worker ---------------- */
  function compute(csvText, fileName) {
    if (worker) { worker.terminate(); worker = null; }
    setProgress(0.01, "Parsing…");
    $("results").style.display = "none";
    try {
      worker = new Worker("worker.js");
    } catch (e) {   // file:// fallback: run on main thread
      try { renderAll(runLocal(csvText), fileName); } catch (err) { setError(err.message); setProgress(null); }
      return;
    }
    worker.onmessage = ev => {
      const m = ev.data;
      if (m.type === "progress") setProgress(m.value, m.label);
      else if (m.type === "error") { setError(m.message); setProgress(null); }
      else if (m.type === "done") { setProgress(null); renderAll(m, fileName); }
    };
    worker.onerror = e => { setError("Worker error: " + e.message); setProgress(null); };
    worker.postMessage({ csvText, opts: opts() });
  }

  function runLocal(csvText) {
    const o = opts();
    const parsed = ReReRe.parseCSV(csvText, o.delimiter);
    const det = ReReRe.detectColumns(parsed.header, parsed.rows);
    if (det.itemCols.length < 6) throw new Error("Fewer than 6 numeric item columns detected.");
    const { M, n, J } = ReReRe.buildMatrix(parsed.rows, det.itemCols);
    const res = ReReRe.ensemble(M, n, J, o);
    const ids = parsed.rows.map((r, i) => det.idCol >= 0 ? r[det.idCol] : "row_" + (i + 1));
    return { type: "done", res, ids, itemNames: det.itemCols.map(c => parsed.header[c]), excluded: det.excluded };
  }

  /* ---------------- rendering ---------------- */
  function renderAll(m, fileName) {
    lastResult = { res: m.res, ids: m.ids, excluded: m.excluded, fileName: fileName || "data.csv" };
    const res = m.res;
    $("results").style.display = "block";

    // diagnostic banner
    const d = res.diagnostic;
    const banner = $("diagBanner");
    banner.className = "banner " + d.level;
    const icon = d.level === "ok" ? ICO.check : ICO.warn;
    const g = res.partnerGate != null ? res.partnerGate : 1;
    const partnerNote = g >= 0.85
      ? " <i>Ensemble partners (LongString, Person-Total) are active.</i>"
      : g > 0.15
        ? " <i>Ensemble partners are partially down-weighted — the item-mean profile is fairly flat, so the score leans on rr.</i>"
        : " <i>Ensemble partners are off — item means are too flat to trust LongString / Person-Total here, so the score is rr alone.</i>";
    banner.innerHTML = "<b>" + icon + " Structure diagnostic: " + d.level.replace("_", " ") +
      "</b> (coherence signal strength " + d.signal_strength.toFixed(2) + "× noise, top-pair mean |r| = " +
      d.top_mean_r.toFixed(2) + ") — " + d.advice + partnerNote;

    // applicability guide (from the ensemble operating-envelope analysis)
    setApplicability(res.params.n, res.params.J);

    // stats (flagging already decided by the engine: automatic or manual)
    $("sN").textContent = res.params.n;
    $("sJ").textContent = res.params.J;
    $("sEst").textContent = Math.round(100 * res.estimatedRate) + "%";
    $("sFlag").textContent = res.nFlagged + " (" + Math.round(100 * res.nFlagged / res.params.n) + "%)";

    updatePrevWarn(res);
    updateParamWarn(res);
    drawScatter(res.p, res.flagged, res.threshold);
    renderTable();
    const ex = m.excluded || [];
    $("exclNote").textContent = ex.length
      ? "Columns set aside (not used as items): " + ex.map(e => e.col + " [" + e.reason + "]").join(", ")
      : "";
  }

  /* Applicability guide — thresholds from the length envelope measured in the
   * controlled simulation, at the SHIPPED cut (severity-weighted MCC, and the
   * ablation that isolates rc's contribution):
   *   >=120 items: wMCC ~.90+, rc's contribution positive with CI excluding 0
   *   >= 60 items: wMCC ~.72, rc's contribution indistinguishable from zero
   *   < 60 items : rc's contribution turns NEGATIVE at 30 items -> partners only.
   * Respondent floors follow the clean-data envelope, whose false-flag tail is
   * worst on small samples. */
  function setApplicability(n, J) {
    const el = $("applicability");
    let cls, msg;
    if (J >= 120 && n >= 200) {
      cls = "ok";
      msg = "<b>" + ICO.check + " Within the validated envelope</b> (" + J + " items, " + n +
        " respondents). This is the regime where the coherence index earns its place: in the " +
        "controlled simulation the ensemble reaches a severity-weighted MCC around 0.90 here, and " +
        "removing <i>rc</i> costs about 0.12 of it.";
    } else if (J >= 60 && n >= 100) {
      cls = "marginal";
      msg = "<b>" + ICO.warn + " Usable, but <i>rc</i> adds little</b> (" + J + " items, " + n +
        " respondents). The ensemble still detects (severity-weighted MCC ≈ 0.72 around 60 items), " +
        "but at this length the coherence index contributes no more than its partners do. " +
        "For its full benefit aim for ≥ 120 items.";
    } else {
      cls = "weak";
      const why = [];
      if (J < 60) why.push("only " + J + " items (≥ 60 to be usable, ≥ 120 recommended)");
      if (n < 100) why.push("only " + n + " respondents (≥ 100 to be usable, ≥ 200 recommended)");
      msg = "<b>" + ICO.warn + " Below the validated envelope:</b> " + why.join(" and ") + ". " +
        "The cross-item structure the ensemble relies on is too thin here — on batteries this short " +
        "the coherence index was measured to <b>subtract</b> from the decision rather than add to it. " +
        "Treat the scores as an <b>exploratory ranking</b>, not a hard flag, and lean on other " +
        "evidence (response times, attention checks).";
    }
    el.className = "banner " + cls;
    el.innerHTML = msg;
  }

  /* Informational note under the results.
   *  - automatic mode, no careless component found -> reassure (nothing flagged).
   *  - manual mode flagging more than the visible low-coherence evidence -> warn
   *    about the top-share cascade (a fixed % is flagged even on clean data). */
  function updatePrevWarn(res) {
    const w = $("prevWarn");
    const N = res.rr.length;
    const evidence = res.rr.filter(v => v <= 1.5).length / N;   // low-coherence rate
    if (res.flagMode === "auto") {
      /* The cut is unconditional, so the question is no longer "did a gate open?"
       * but "is this flagged share more than clean data of this size produces by
       * chance?" — which is what the measured envelope answers. */
      const env = res.cleanEnvelope;
      const pct = x => (100 * x).toFixed(x < 0.1 ? 1 : 0) + "%";
      if (res.flagStatus === "none-past-cut") {
        w.className = "banner ok"; w.style.display = "";
        w.innerHTML = ICO.check + " <b>Nobody reaches the cutoff.</b> No respondent scores more than " +
          "2.5&sigma; above the attentive cluster, so nothing is flagged — the expected result on clean " +
          "(or already-cleaned) data. Carelessness that shifts a profile slightly, without producing " +
          "extreme responses, can also look like this: the continuous ranking below still orders " +
          "respondents, and an expected rate can be set under Advanced settings.";
      } else if (env && res.estimatedRate <= env.p95) {
        w.className = "banner marginal"; w.style.display = "";
        w.innerHTML = ICO.warn + " <b>This many flags is within what clean data of this size produces by " +
          "chance.</b> On simulated <i>careless-free</i> samples with about " + res.params.J + " items and " +
          res.params.n + " respondents, the same rule flags " + pct(env.mean) + " on average and up to " +
          pct(env.p95) + " in the worst 1-in-20 — against the " + pct(res.estimatedRate) + " flagged here. " +
          "Treat these flags as unconfirmed: inspect them individually rather than removing them in bulk. " +
          "The envelope narrows quickly with more items and more respondents.";
      } else { w.style.display = "none"; }
      return;
    }
    // manual override
    const rate = res.params.prevalence || 0;
    if (rate > evidence * 1.5 + 0.02) {
      const pct = Math.round(rate * 100), ev = Math.round(evidence * 100);
      w.className = "banner weak"; w.style.display = "";
      w.innerHTML = ICO.warn + " <b>Some of these flags may be false positives.</b> You set a manual rate of " +
        pct + "%, but only about <b>" + ev + "%</b> of respondents show clearly low response coherence. A " +
        "top-share rule marks that share even on clean data, so if this file was already cleaned don't treat " +
        "the extras as new careless. Leave the rate blank for automatic calibration, or inspect borderline rows.";
    } else { w.style.display = "none"; }
  }

  /* Out-of-envelope PARAMETER warning: fires when the run used settings outside the
   * ranges we validated (a custom sensitivity z, an extreme corProp, too few
   * permutations, etc.). Freedom is allowed, but the user is told the numbers no
   * longer carry the study's guarantees. */
  function updateParamWarn(res) {
    const w = $("paramWarn"); if (!w) return;
    const p = res.params, issues = [];
    const s = p.sensitivity;
    if (typeof s === "number" && s !== 2.5)
      issues.push("a custom sensitivity cut of <b>z = " + s + "</b> instead of the shipped z = 2.5"
        + ((s < 1.5 || s > 2.5) ? " (outside the 1.5–2.5 range we evaluated at all)" : ""));
    if (p.corProp != null && (p.corProp < 0.02 || p.corProp > 0.20))
      issues.push("<b>corProp = " + p.corProp + "</b> (studied range 0.03–0.20)");
    if (p.iterations != null && p.iterations < 100)
      issues.push("only <b>" + p.iterations + "</b> permutations (≥ 100 for a stable baseline)");
    if (p.minPairs != null && p.minPairs < 10)
      issues.push("<b>minPairs = " + p.minPairs + "</b> (≥ 15 recommended)");
    if (!issues.length) { w.style.display = "none"; return; }
    w.className = "banner weak"; w.style.display = "";
    w.innerHTML = ICO.warn + " <b>Running outside the studied limits.</b> This analysis used " +
      issues.join("; ") + ". These settings were not covered by our validation, so the detection " +
      "quality we report may not hold here — treat the output as <b>exploratory</b> and revert to the " +
      "defaults for validated behaviour.";
  }

  let sortKey = "p", sortAsc = false;
  function renderTable() {
    const { res, ids } = lastResult;
    const rows = ids.map((id, i) => ({
      id, p: res.p[i],
      p_se: res.etaSe ? res.etaSe[i] * res.p[i] * (1 - res.p[i]) : 0,
      borderline: res.borderline ? res.borderline[i] : false,
      flagged: res.flagged[i], rr: res.rr[i], longstring: res.longstring[i],
      person_total: res.person_total[i], irv: res.irv[i], d2: res.d2[i], n_missing: res.nMissing[i]
    }));
    rows.sort((a, b) => {
      const va = a[sortKey], vb = b[sortKey];
      const c = (typeof va === "string") ? String(va).localeCompare(String(vb)) : (va - vb);
      return sortAsc ? c : -c;
    });
    const cols = [["id", "ID"], ["p", "P(careless)"], ["flagged", "flagged"], ["rr", "rr"],
                  ["longstring", "LongString"], ["person_total", "Person-Total"],
                  ["irv", "IRV"], ["d2", "D²"], ["n_missing", "missing"]];
    const cap = 500;
    let html = "<thead><tr>" + cols.map(c =>
      "<th data-k='" + c[0] + "'>" + c[1] + (sortKey === c[0] ? (sortAsc ? " ▲" : " ▼") : "") + "</th>"
    ).join("") + "</tr></thead><tbody>";
    rows.slice(0, cap).forEach(r => {
      html += "<tr" + (r.flagged ? " class='flag'" : "") + "><td>" + escapeHtml(String(r.id)) + "</td>" +
        "<td>" + r.p.toFixed(3) + " <span class='se'>±" + r.p_se.toFixed(3) + "</span></td>" +
        "<td class='flagcell'>" + (r.flagged ? ICO.flag + " careless" + (r.borderline ? " <span class='se'>(borderline)</span>" : "") : "") + "</td>" +
        "<td>" + r.rr.toFixed(2) + "</td><td>" + r.longstring + "</td>" +
        "<td>" + r.person_total.toFixed(2) + "</td>" +
        "<td>" + r.irv.toFixed(3) + "</td><td>" + r.d2.toFixed(1) + "</td><td>" + r.n_missing + "</td></tr>";
    });
    html += "</tbody>";
    const tbl = $("tbl");
    tbl.innerHTML = html;
    tbl.querySelectorAll("th").forEach(th => th.addEventListener("click", () => {
      const k = th.getAttribute("data-k");
      if (sortKey === k) sortAsc = !sortAsc; else { sortKey = k; sortAsc = true; }
      renderTable();
    }));
    $("tblNote").textContent = rows.length > cap
      ? "Showing first " + cap + " of " + rows.length + " rows (sorted). The download contains all rows."
      : "";
  }

  function escapeHtml(s) {
    return s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;").replace(/"/g, "&quot;");
  }

  /* Per-respondent dot plot: every respondent is a point, ranked left→right by
   * careless probability (y). Flagged points are red, the rest slate-blue; a
   * dashed line marks the cutoff. Shows the separation (or lack of it) directly. */
  function drawScatter(p, flagged, thr) {
    const cv = $("hist");
    const W = cv.clientWidth || 900, H = 230;
    cv.width = W * devicePixelRatio; cv.height = H * devicePixelRatio;
    const g = cv.getContext("2d");
    g.scale(devicePixelRatio, devicePixelRatio);
    g.clearRect(0, 0, W, H);
    const n = p.length; if (!n) return;
    const idx = [...Array(n).keys()].filter(i => Number.isFinite(p[i])).sort((a, b) => p[a] - p[b]);
    const padL = 40, padR = 14, padB = 28, padT = 12;
    const plotW = W - padL - padR, plotH = H - padT - padB;
    const X = r => padL + plotW * (n > 1 ? r / (n - 1) : 0.5);
    const Y = v => padT + plotH * (1 - v);                 // p in [0,1]
    // horizontal gridlines + y labels
    g.font = '11px "Iowan Old Style", Palatino, Georgia, serif';
    [0, 0.25, 0.5, 0.75, 1].forEach(t => {
      const y = Y(t);
      g.strokeStyle = "#ece9e1"; g.lineWidth = 1;
      g.beginPath(); g.moveTo(padL, y); g.lineTo(W - padR, y); g.stroke();
      g.fillStyle = "#5f6b7a"; g.fillText(t.toFixed(2), 8, y + 4);
    });
    // cutoff line
    if (Number.isFinite(thr) && thr >= 0 && thr <= 1) {
      const y = Y(thr);
      g.strokeStyle = "#8f3535"; g.lineWidth = 1.5; g.setLineDash([5, 4]);
      g.beginPath(); g.moveTo(padL, y); g.lineTo(W - padR, y); g.stroke(); g.setLineDash([]);
      g.fillStyle = "#8f3535"; g.fillText("cutoff " + thr.toFixed(2), W - padR - 96, y - 5);
    }
    // points
    for (let r = 0; r < idx.length; r++) {
      const i = idx[r];
      g.globalAlpha = 0.82;
      g.fillStyle = flagged[i] ? "#b0553f" : "#5f7c9e";
      g.beginPath(); g.arc(X(r), Y(p[i]), 2.6, 0, 2 * Math.PI); g.fill();
    }
    g.globalAlpha = 1;
    g.fillStyle = "#5f6b7a"; g.font = 'italic 11px "Iowan Old Style", Palatino, Georgia, serif';
    g.fillText("respondents ranked by careless probability →", padL, H - 8);
  }

  /* ---------------- downloads ---------------- */
  function resultsCsv(onlyFlagged) {
    const { res, ids } = lastResult;
    const header = "id,careless_prob,careless_prob_se,flagged,borderline,rc,longstring,person_total,irv,d2,n_missing";
    const lines = [header];
    ids.forEach((id, i) => {
      if (onlyFlagged && !res.flagged[i]) return;
      const idq = /[",;\n]/.test(String(id)) ? '"' + String(id).replace(/"/g, '""') + '"' : id;
      const pse = res.etaSe ? res.etaSe[i] * res.p[i] * (1 - res.p[i]) : 0;
      lines.push([idq, res.p[i].toFixed(4), pse.toFixed(4), res.flagged[i] ? 1 : 0,
                  (res.borderline && res.borderline[i]) ? 1 : 0, res.rr[i].toFixed(4),
                  res.longstring[i], res.person_total[i].toFixed(4), res.irv[i].toFixed(4),
                  res.d2[i].toFixed(3), res.nMissing[i]].join(","));
    });
    return lines.join("\n");
  }
  function download(name, text) {
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([text], { type: "text/csv" }));
    a.download = name;
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 3000);
  }
  $("dlBtn").addEventListener("click", () => {
    if (!lastResult) return;
    download(lastResult.fileName.replace(/\.\w+$/, "") + "_rerere.csv", resultsCsv(false));
  });
  $("dlFlagBtn").addEventListener("click", () => {
    if (!lastResult) return;
    download(lastResult.fileName.replace(/\.\w+$/, "") + "_flagged.csv", resultsCsv(true));
  });

  /* live re-flag when the manual rate changes — no recompute, only the cut moves.
   * Blank field => automatic attentive-mode calibration; a value => top-share. */
  function reflag(res) {
    const mr = manualRate();
    if (mr != null) {
      const thr = applyPrevalence(res, mr);
      res.flagMode = "manual";
      res.flagStatus = "manual";
      res.estimatedRate = res.nFlagged / res.params.n;
    } else {
      const af = ReReRe.autoFlag(res.eta, sensSetting());
      res.flagged = af.flagged;
      res.nFlagged = af.flagged.filter(Boolean).length;
      res.threshold = af.twoComp ? 1 / (1 + Math.exp(-af.etaThreshold)) : Infinity;
      res.estimatedRate = af.estRate;
      res.flagMode = "auto";
      res.twoComponent = af.twoComp;
      res.flagStatus = af.status;
      res.params.prevalence = null;
      res.params.sensitivity = af.sensitivity;
    }
    $("sEst").textContent = Math.round(100 * res.estimatedRate) + "%";
    $("sFlag").textContent = res.nFlagged + " (" + Math.round(100 * res.nFlagged / res.params.n) + "%)";
    updatePrevWarn(res);
    updateParamWarn(res);
    drawScatter(res.p, res.flagged, res.threshold);
    renderTable();
  }
  $("pPrev").addEventListener("input", () => { if (lastResult) reflag(lastResult.res); });
  $("pSensZ").addEventListener("input", () => { if (lastResult) reflag(lastResult.res); });
})();
