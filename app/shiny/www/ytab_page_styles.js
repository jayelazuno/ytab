// Bootstrap 3 is loaded for the whole app through shiny::bootstrapLib() in the head: Shiny's dialogs, the shinyFiles
// directory chooser and every Shiny input need its stylesheet and its JavaScript wherever they appear.
//
// Bootstrap sets `html { font-size: 10px }`, and YTAB's own rules are written in rem. The workspace is designed
// against that 10px root; the landing page is designed against the browser's usual 16px one. Marking the document
// while the landing page is shown lets ytab.css restore the 16px root for it alone, so each page keeps the
// proportions it was designed with.
(function () {
  function update() {
    var landing = document.querySelector(".ytab-landing") !== null;
    document.documentElement.classList.toggle("ytab-landing-root", landing);
  }
  // The Maximize/Restore button in a dialog's title row. The class goes on the dialog, so the size lasts as long as
  // the dialog does and no server round trip is needed; a dialog opened again starts at its usual size.
  document.addEventListener("click", function (event) {
    var target = event.target;
    var button = target && target.closest ? target.closest("[data-ytab-modal-size]") : null;
    if (!button) return;
    var dialog = button.closest(".modal-dialog");
    if (!dialog) return;
    var maximized = dialog.classList.toggle("ytab-modal-max");
    button.textContent = maximized ? "Restore" : "Maximize";
    button.setAttribute("aria-label", maximized ? "Return the dialog to its usual size" : "Fill the window");
  });

  new MutationObserver(update).observe(document.documentElement, { childList: true, subtree: true });
  document.addEventListener("DOMContentLoaded", update);
  update();
})();
