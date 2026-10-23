(function () {
  "use strict";

  const menus = {
    project: {
      nested_id: "project_tabs",
      items: [
        ["Samples", "samples"],
        ["Logs & Provenance", "provenance"],
        ["Help", "help"]
      ]
    },
    run: {
      nested_id: "run_tabs",
      items: [
        ["Samples", "samples"],
        ["Reference", "reference"],
        ["Mapping", "mapping"],
        ["Run Profiles", "profiles"]
      ]
    },
    quality_control: {
      nested_id: "qc_tabs",
      items: [
        ["Overview", "overview"],
        ["Reads & Mapping", "reads_mapping"],
        ["Insertions & Genes", "insertions"],
        ["Saturation & Bias", "saturation_bias"],
        ["Library Agreement", "agreement"]
      ]
    },
    normalization: {
      nested_id: "normalization_tabs",
      items: [
        ["Site Normalization", "site"]
      ]
    },
    essentiality: {
      nested_id: "essentiality_tabs",
      items: [
        ["Run", "run"],
        ["Results", "results"]
      ]
    },
    fitness: {
      nested_id: "fitness_tabs",
      items: [
        ["Design", "fitness_design"],
        ["Treated vs Control: Run", "fitness_run"],
        ["Treated vs Control: Results", "fitness_results"],
        ["Multi-condition", "multicondition"],
        ["Genetic Interactions", "gi"]
      ]
    },
    gene_explorer: {
      nested_id: "gene_explorer_tabs",
      items: [
        ["Genome Browser", "browser"],
        ["Gene & Domain Insertions", "gene_domain"],
        ["GO Enrichment", "go"],
        ["Across Species", "comparative"]
      ]
    },
    screen_coverage: {
      nested_id: "screen_coverage_tabs",
      items: [["Coverage Planner", "planner"], ["Mutant Identification", "mutants"]]
    },
  };

  const nestedByTop = {};
  Object.keys(menus).forEach(function (top) {
    nestedByTop[top] = menus[top].nested_id;
  });

  let applyingServerState = false;

  function closeMenus(except) {
    document.querySelectorAll(".ytab-nav-menu-host.open").forEach(function (host) {
      if (host !== except) {
        host.classList.remove("open");
        const toggle = host.querySelector(".ytab-nav-menu-caret");
        if (toggle) toggle.setAttribute("aria-expanded", "false");
      }
    });
  }

  function clearHoverSuppression(host) {
    if (host) delete host.dataset.ytabSuppressHover;
  }

  function suppressHoverUntilLeave(host) {
    if (!host) return;
    host.dataset.ytabSuppressHover = "true";
    setOpen(host, false);
  }

  function setOpen(host, open) {
    closeMenus(open ? host : null);
    host.classList.toggle("open", open);
    const toggle = host.querySelector(".ytab-nav-menu-caret");
    if (toggle) toggle.setAttribute("aria-expanded", open ? "true" : "false");
  }

  function sendShinyValue(inputId, payload) {
    if (!window.Shiny || typeof window.Shiny.setInputValue !== "function") return;
    window.Shiny.setInputValue(inputId, JSON.stringify(payload), { priority: "event" });
  }

  function sendMenuTarget(top, nestedId, nestedValue) {
    sendShinyValue("ytab_nav_menu", {
      top: top,
      nested_id: nestedId,
      nested_value: nestedValue,
      nonce: Date.now()
    });
  }

  function parseHash() {
    const raw = window.location.hash ? window.location.hash.replace(/^#/, "") : "";
    if (!raw) return null;
    const params = new URLSearchParams(raw);
    const top = params.get("tab") || params.get("top") || "";
    if (!top) return null;
    const nestedValue = params.get("subtab") || params.get("nested_value") || "";
    return {
      top: top,
      nested_id: params.get("nested_id") || nestedByTop[top] || "",
      nested_value: nestedValue,
      results_view: params.get("results_view") || "",
      project_id: params.get("project") || params.get("project_id") || "",
      nonce: Date.now()
    };
  }

  function hashForState(state) {
    if (!state || !state.top) return "";
    const params = new URLSearchParams();
    if (state.project_id) params.set("project", state.project_id);
    params.set("tab", state.top);
    if (state.nested_id) params.set("nested_id", state.nested_id);
    if (state.nested_value) params.set("subtab", state.nested_value);
    if (state.results_view) params.set("results_view", state.results_view);
    return "#" + params.toString();
  }

  function writeHash(state) {
    const next = hashForState(state);
    if (!next || next === window.location.hash) return;
    applyingServerState = true;
    if (state && state.replace) {
      window.history.replaceState(null, "", next);
    } else {
      window.history.pushState(null, "", next);
    }
    window.setTimeout(function () {
      applyingServerState = false;
    }, 0);
  }

  function sendHashLocation() {
    if (applyingServerState) return;
    const target = parseHash();
    if (!target) return;
    sendShinyValue("ytab_nav_hash", target);
  }

  function restoreFromBrowserHistory() {
    if (applyingServerState) return;
    const target = parseHash();
    if (!target) return;
    sendShinyValue("ytab_nav_hash", target);
    window.setTimeout(function () {
      clientNavigate(target.top, target.nested_value);
    }, 180);
  }

  function clickTab(value) {
    const candidates = Array.from(document.querySelectorAll('a[data-value="' + value + '"], button[data-value="' + value + '"]'));
    const visible = candidates.find(function (node) {
      return node.offsetParent !== null;
    }) || candidates[0];
    if (visible) visible.click();
    return Boolean(visible);
  }

  function clientNavigate(top, nestedValue) {
    clickTab(top);
    window.setTimeout(function () {
      if (nestedValue) clickTab(nestedValue);
    }, 80);
  }

  function installMenus() {
    Object.keys(menus).forEach(function (top) {
      const config = menus[top];
      const link = document.querySelector('.navbar-nav > li > a[data-value="' + top + '"]');
      if (!link) return;
      const host = link.closest("li");
      if (!host || host.dataset.ytabMenuReady === "true") return;

      host.dataset.ytabMenuReady = "true";
      // Not "dropdown": Shiny's tab input reads the active tab with li:not(.dropdown).active, so a top tab marked
      // as a dropdown would report no value and the server would never see which module is open.
      host.classList.add("ytab-nav-menu-host");

      const caret = document.createElement("button");
      caret.type = "button";
      caret.className = "ytab-nav-menu-caret";
      caret.setAttribute("aria-label", "Open " + link.textContent.trim() + " menu");
      caret.setAttribute("aria-expanded", "false");
      caret.textContent = "▾";

      const menu = document.createElement("ul");
      menu.className = "dropdown-menu ytab-nav-submenu";

      config.items.forEach(function (entry) {
        const li = document.createElement("li");
        const item = document.createElement("a");
        item.href = "#";
        item.textContent = entry[0];
        item.addEventListener("click", function (event) {
          event.preventDefault();
          event.stopPropagation();
          clientNavigate(top, entry[1]);
          sendMenuTarget(top, config.nested_id, entry[1]);
          suppressHoverUntilLeave(host);
          item.blur();
          if (document.activeElement && typeof document.activeElement.blur === "function") {
            document.activeElement.blur();
          }
        });
        li.appendChild(item);
        menu.appendChild(li);
      });

      host.addEventListener("mouseenter", function () {
        if (host.dataset.ytabSuppressHover === "true") return;
        setOpen(host, true);
      });

      host.addEventListener("mouseleave", function () {
        clearHoverSuppression(host);
        setOpen(host, false);
      });

      link.addEventListener("click", function () {
        suppressHoverUntilLeave(host);
        if (typeof link.blur === "function") link.blur();
      });

      caret.addEventListener("click", function (event) {
        event.preventDefault();
        event.stopPropagation();
        clearHoverSuppression(host);
        setOpen(host, !host.classList.contains("open"));
      });

      link.addEventListener("keydown", function (event) {
        if (event.key === "ArrowDown") {
          event.preventDefault();
          setOpen(host, true);
          const first = menu.querySelector("a");
          if (first) first.focus();
        }
      });

      menu.addEventListener("keydown", function (event) {
        const links = Array.from(menu.querySelectorAll("a"));
        const index = links.indexOf(document.activeElement);
        if (event.key === "Escape") {
          event.preventDefault();
          setOpen(host, false);
          caret.focus();
        } else if (event.key === "ArrowDown") {
          event.preventDefault();
          links[(index + 1 + links.length) % links.length].focus();
        } else if (event.key === "ArrowUp") {
          event.preventDefault();
          links[(index - 1 + links.length) % links.length].focus();
        }
      });

      host.appendChild(caret);
      host.appendChild(menu);
    });
  }

  document.addEventListener("click", function (event) {
    if (!event.target.closest(".ytab-nav-menu-host")) closeMenus();
  });

  window.addEventListener("hashchange", restoreFromBrowserHistory);
  window.addEventListener("popstate", restoreFromBrowserHistory);
  document.addEventListener("shiny:connected", sendHashLocation);
  if (window.Shiny && typeof window.Shiny.addCustomMessageHandler === "function") {
    window.Shiny.addCustomMessageHandler("ytab-nav-state", writeHash);
  } else {
    document.addEventListener("shiny:connected", function () {
      window.Shiny.addCustomMessageHandler("ytab-nav-state", writeHash);
    }, { once: true });
  }

  document.addEventListener("DOMContentLoaded", installMenus);
  new MutationObserver(installMenus).observe(document.documentElement, {
    childList: true,
    subtree: true
  });
})();
