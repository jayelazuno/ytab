/* Saving the figure that is on screen.
 *
 * The download buttons used to be served by a second, server-side ggplot rebuild of each plot. Keeping two
 * builders agreeing proved impossible in practice: panels, scales, chromosome names and whole layers drifted
 * apart without anything failing. Plotly already holds the figure the user is looking at, and Plotly.toImage
 * renders exactly that -- the same call behind the mode bar's camera icon -- so the button captures the real
 * figure instead of approximating it.
 *
 * PNG is rasterised at `scale` (3 is roughly 300 dpi for a 12-inch figure). PDF is not something plotly.js can
 * write, so vector output is SVG, which is what a manuscript wants anyway.
 */
(function () {
  "use strict";

  function graphOf(outputId) {
    var host = document.getElementById(outputId);
    if (!host) return null;
    // plotlyOutput renders into a child div; the element carrying _fullLayout is the figure itself.
    if (host._fullLayout) return host;
    return host.querySelector(".js-plotly-plot");
  }

  function capture(message) {
    var graph = graphOf(message.target);
    if (!graph || !window.Plotly) {
      if (window.Shiny) Shiny.setInputValue(message.status, "missing", { priority: "event" });
      return;
    }
    var box = graph.getBoundingClientRect();
    Plotly.toImage(graph, {
      format: message.format === "svg" ? "svg" : "png",
      width: Math.max(900, Math.round(box.width)),
      height: Math.max(400, Math.round(box.height)),
      scale: message.format === "svg" ? 1 : (message.scale || 3)
    }).then(function (url) {
      var link = document.createElement("a");
      link.href = url;
      link.download = message.filename;
      document.body.appendChild(link);
      link.click();
      document.body.removeChild(link);
      if (window.Shiny) Shiny.setInputValue(message.status, "saved", { priority: "event" });
    }).catch(function () {
      if (window.Shiny) Shiny.setInputValue(message.status, "failed", { priority: "event" });
    });
  }

  function register() {
    if (!window.Shiny || typeof window.Shiny.addCustomMessageHandler !== "function") return false;
    window.Shiny.addCustomMessageHandler("ytab-capture-plot", capture);
    return true;
  }

  if (!register()) document.addEventListener("DOMContentLoaded", register);
})();
