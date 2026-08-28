(function () {
  let browser = null;
  let lastConfig = null;
  let lastLocus = "";
  let loadSerial = 0;

  function warning(message) {
    const node = document.getElementById("ytab_igv_warning");
    if (!node) return;
    if (!message) {
      node.style.display = "none";
      node.textContent = "";
      return;
    }
    node.style.display = "block";
    node.textContent = message;
  }

  function container() {
    return document.getElementById("ytab_igv_browser");
  }

  function genomeBrowserIsVisible() {
    const igvDiv = container();
    if (!igvDiv) return false;
    const rect = igvDiv.getBoundingClientRect();
    return rect.width > 0 && rect.height > 0;
  }

  function normalizeTrack(track) {
    const cleaned = Object.assign({}, track);
    if (cleaned.url && !/^https?:\/\//i.test(cleaned.url) && !cleaned.url.startsWith("/")) {
      cleaned.url = "/" + cleaned.url;
    }
    return cleaned;
  }

  function normalizeConfig(config) {
    const next = Object.assign({}, config || {});
    next.genome = Object.assign({}, next.genome || {});
    ["fastaURL", "indexURL", "aliasURL", "cytobandURL", "chromSizesURL"].forEach(function (key) {
      if (next.genome[key] && !/^https?:\/\//i.test(next.genome[key]) && !next.genome[key].startsWith("/")) {
        next.genome[key] = "/" + next.genome[key];
      }
    });
    next.genome.tracks = (next.genome.tracks || []).map(normalizeTrack);
    next.tracks = (next.tracks || []).map(normalizeTrack);
    return next;
  }

  function currentLocusString() {
    if (browser && browser.currentLoci) {
      const loci = browser.currentLoci();
      if (typeof loci === "string" && loci) return loci;
      if (loci && loci.length) {
        if (typeof loci[0] === "string") return loci[0];
        const chr = loci[0].chr || loci[0].chrName || loci[0].chromosome || "";
        const start = Math.max(1, Math.round(loci[0].start || 1));
        const end = Math.max(start + 1, Math.round(loci[0].end || start + 1));
        if (chr) return chr + ":" + start + "-" + end;
      }
    }
    if (lastLocus) return lastLocus;
    const igvDiv = container();
    if (igvDiv) {
      const inputs = Array.from(igvDiv.querySelectorAll("input"));
      const hit = inputs
        .map(function (input) { return String(input.value || ""); })
        .find(function (value) {
          return /^([^:]+):[0-9,]+-[0-9,]+$/.test(value.trim());
        });
      if (hit) return hit;
    }
    const frames = browser && (browser.referenceFrameList || browser.referenceFrames);
    const frame = frames && frames.length ? frames[0] : browser && browser.referenceFrame;
    if (frame) {
      const chr = frame.chrName || frame.chr || frame.chromosome || "";
      const start = Math.max(1, Math.round((frame.start || 0) + 1));
      const bpPerPixel = frame.bpPerPixel || 1;
      const width = (container() && container().getBoundingClientRect().width) || 1000;
      const end = Math.max(start + 1, Math.round(start + bpPerPixel * width));
      if (chr) return chr + ":" + start + "-" + end;
    }
    return "";
  }

  function pan(direction) {
    const locus = currentLocusString();
    const match = String(locus).replace(/,/g, "").match(/^([^:]+):([0-9]+)-([0-9]+)$/);
    if (!match || !browser) return;
    const chr = match[1];
    const start = parseInt(match[2], 10);
    const end = parseInt(match[3], 10);
    const width = Math.max(1, end - start + 1);
    const delta = Math.max(1, Math.round(width * 0.35)) * direction;
    const nextStart = Math.max(1, start + delta);
    const nextEnd = Math.max(nextStart + 1, end + delta);
    const nextLocus = chr + ":" + nextStart + "-" + nextEnd;
    if (typeof browser.search === "function") {
      try {
        lastLocus = nextLocus;
        const searched = browser.search(nextLocus);
        if (searched && typeof searched.catch === "function") {
          searched.catch(function () {
            if (typeof browser.goto === "function") browser.goto(chr, nextStart - 1, nextEnd);
          });
        }
        return;
      } catch (ignored) {}
    }
    if (typeof browser.goto === "function") {
      lastLocus = nextLocus;
      browser.goto(chr, nextStart - 1, nextEnd);
    }
  }

  function locusFromReferenceFrames(referenceFrameList) {
    const frames = Array.isArray(referenceFrameList) ? referenceFrameList : [referenceFrameList];
    const loci = frames
      .filter(function (frame) { return frame; })
      .map(function (frame) {
        if (typeof frame.getLocusString === "function") return frame.getLocusString();
        const chr = frame.chrName || frame.chr || frame.chromosome || "";
        const start = Math.max(1, Math.round((frame.start || 0) + 1));
        const bpPerPixel = frame.bpPerPixel || 1;
        const width = (container() && container().getBoundingClientRect().width) || 1000;
        const end = Math.max(start + 1, Math.round(start + bpPerPixel * width));
        return chr ? chr + ":" + start + "-" + end : "";
      })
      .filter(function (locus) { return locus; });
    return loci.length ? loci.join(" ") : "";
  }

  function attachLocusChangeHandler() {
    if (!browser || typeof browser.on !== "function" || browser.__ytabLocusChangeAttached) return;
    browser.__ytabLocusChangeAttached = true;
    browser.on("locuschange", function (referenceFrameList) {
      const locus = locusFromReferenceFrames(referenceFrameList);
      if (locus) lastLocus = locus;
    });
  }

  function handlePanKeydown(event) {
    if (!genomeBrowserIsVisible() || !browser) return;
    if (event.key !== "ArrowLeft" && event.key !== "ArrowRight") return;
    const igvDiv = container();
    const active = document.activeElement;
    const target = event.target || active;
    const tag = active && active.tagName ? active.tagName.toLowerCase() : "";
    const activeInsideBrowser = igvDiv && active && igvDiv.contains(active);
    const targetInsideBrowser = igvDiv && target && igvDiv.contains(target);
    const editingOutsideBrowser = !activeInsideBrowser &&
      (["input", "textarea", "select"].includes(tag) || active && active.isContentEditable);
    if (editingOutsideBrowser) return;
    if (!activeInsideBrowser && !targetInsideBrowser && !document.body.classList.contains("ytab-genome-browser-active")) return;
    event.preventDefault();
    if (typeof event.stopImmediatePropagation === "function") event.stopImmediatePropagation();
    else event.stopPropagation();
    pan(event.key === "ArrowLeft" ? -1 : 1);
  }

  function attachKeyboardPan() {
    const igvDiv = container();
    if (!igvDiv || igvDiv.dataset.ytabKeyboardPan === "1") return;
    igvDiv.dataset.ytabKeyboardPan = "1";
    igvDiv.tabIndex = 0;
    igvDiv.addEventListener("click", function () { igvDiv.focus({ preventScroll: true }); });
    igvDiv.addEventListener("keydown", handlePanKeydown, true);
  }

  function attachDocumentKeyboardPan() {
    if (document.body.dataset.ytabIgvKeyboardPan === "1") return;
    document.body.dataset.ytabIgvKeyboardPan = "1";
    document.addEventListener("keydown", handlePanKeydown, true);
    window.addEventListener("keydown", handlePanKeydown, true);
  }

  function attachPanButtons() {
    if (document.body.dataset.ytabIgvPanButtons === "1") return;
    document.body.dataset.ytabIgvPanButtons = "1";
    document.addEventListener("click", function (event) {
      const target = event.target;
      const left = target && target.closest ? target.closest("#genome_browser_pan_left") : null;
      const right = target && target.closest ? target.closest("#genome_browser_pan_right") : null;
      if (!left && !right) return;
      event.preventDefault();
      pan(left ? -1 : 1);
    }, true);
  }

  async function createOrReload(config) {
    const serial = ++loadSerial;
    const igvDiv = container();
    if (!igvDiv) return;
    document.body.classList.add("ytab-genome-browser-active");
    attachKeyboardPan();
    attachDocumentKeyboardPan();
    attachPanButtons();
    if (typeof igv === "undefined" || !igv.createBrowser) {
      warning("IGV.js is not available. Check network access to the IGV JavaScript dependency.");
      return;
    }
    warning("");
    lastConfig = normalizeConfig(config);
    lastLocus = lastConfig.locus || "";
    const igvConfig = {
      reference: lastConfig.genome,
      locus: lastConfig.locus || undefined,
      tracks: lastConfig.tracks || []
    };
    if (browser) {
      try {
        await browser.loadSessionObject(igvConfig);
        attachLocusChangeHandler();
        return;
      } catch (e) {
        try { igv.removeBrowser(browser); } catch (ignored) {}
        browser = null;
      }
    }
    try { igv.removeAllBrowsers && igv.removeAllBrowsers(); } catch (ignored) {}
    igvDiv.innerHTML = "";
    const nextBrowser = await igv.createBrowser(igvDiv, igvConfig);
    if (serial !== loadSerial) {
      try { igv.removeBrowser && igv.removeBrowser(nextBrowser); } catch (ignored) {}
      return;
    }
    browser = nextBrowser;
    attachLocusChangeHandler();
    attachKeyboardPan();
    attachPanButtons();
  }

  async function setTracks(tracks) {
    if (!browser) {
      if (lastConfig) {
        lastConfig.tracks = tracks || [];
        await createOrReload(lastConfig);
      }
      return;
    }
    const normalized = (tracks || []).map(normalizeTrack);
    const currentTracks = (browser.trackViews || [])
      .map(function (view) { return view.track; })
      .filter(function (track) { return track && track.removable !== false; });
    for (const track of currentTracks) {
      try { browser.removeTrack(track); } catch (e) {}
    }
    for (const track of normalized) {
      await browser.loadTrack(track);
    }
  }

  function gotoLocus(locus) {
    if (!browser || !locus) return;
    lastLocus = String(locus);
    browser.search(String(locus));
  }

  function attachHandlers() {
    if (!window.Shiny || !Shiny.addCustomMessageHandler) {
      setTimeout(attachHandlers, 100);
      return;
    }
    Shiny.addCustomMessageHandler("ytab_igv_init", createOrReload);
    Shiny.addCustomMessageHandler("ytab_igv_set_tracks", setTracks);
    Shiny.addCustomMessageHandler("ytab_igv_goto", gotoLocus);
    Shiny.addCustomMessageHandler("ytab_igv_pan", function (direction) {
      pan(Number(direction) < 0 ? -1 : 1);
    });
    Shiny.addCustomMessageHandler("ytab_igv_warning", warning);
  }

  attachPanButtons();
  attachHandlers();
})();
