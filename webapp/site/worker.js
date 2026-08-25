/* worker.js — runs the ReReRe ensemble computation off the main thread. */
importScripts("reread.js");

onmessage = function (ev) {
  const { csvText, opts } = ev.data;
  try {
    postMessage({ type: "progress", value: 0.02, label: "Parsing CSV…" });
    const parsed = ReReRe.parseCSV(csvText, opts.delimiter || undefined);
    const det = ReReRe.detectColumns(parsed.header, parsed.rows);
    if (det.itemCols.length < 6)
      throw new Error("Fewer than 6 numeric item columns detected — check the delimiter (Advanced settings) and that items are numeric.");
    const { M, n, J } = ReReRe.buildMatrix(parsed.rows, det.itemCols);
    postMessage({ type: "progress", value: 0.05, label: "Computing correlations (" + n + " × " + J + ")…" });
    const res = ReReRe.ensemble(M, n, J, Object.assign({}, opts, {
      onProgress: p => postMessage({
        type: "progress", value: p,
        label: p < 0.16 ? "Computing correlations…" : "Permutation baseline… " + Math.round(p * 100) + "%"
      })
    }));
    const ids = parsed.rows.map((r, i) => det.idCol >= 0 ? r[det.idCol] : "row_" + (i + 1));
    postMessage({ type: "done", res, ids,
                  itemNames: det.itemCols.map(c => parsed.header[c]),
                  excluded: det.excluded });
  } catch (e) {
    postMessage({ type: "error", message: e.message || String(e) });
  }
};
