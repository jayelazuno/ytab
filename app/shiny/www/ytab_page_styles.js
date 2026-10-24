// The workspace (navbarPage) loads Bootstrap 3, and its stylesheets stay in the page after "Switch project".
// The landing page is designed without Bootstrap, so switch those stylesheets off whenever the landing page is shown
// and back on for the workspace; the landing page then looks the same as when the app first opens. The media attribute
// is used rather than `disabled`: a sheet disabled while it is still loading is not applied when re-enabled.
(function () {
  function bootstrapSheets() {
    return Array.prototype.filter.call(document.querySelectorAll('link[rel="stylesheet"]'), function (link) {
      return /bootstrap-[0-9.]+\/(css|accessibility)\//.test(link.getAttribute("href") || "");
    });
  }
  function update() {
    var landing = document.querySelector(".ytab-landing") !== null;
    var media = landing ? "not all" : "all";
    bootstrapSheets().forEach(function (link) { if (link.media !== media) link.media = media; });
  }
  new MutationObserver(update).observe(document.documentElement, { childList: true, subtree: true });
  document.addEventListener("DOMContentLoaded", update);
})();
