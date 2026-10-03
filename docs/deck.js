// deck.js: the slide deck.
// Scales every 1600 x 900 slide to fit the window next to the rail, wires the tabs, marks the current page in the
// rail, and tells the page's own code when the current page or the scale changes ("framechange", "slidescale").
(() => {
  const W = 1600, H = 900;
  const root = document.documentElement;
  const rail = document.querySelector(".rail");

  // ---- fit: one scale for every slide, the largest that fits beside the rail ----
  let scale = 0;
  function fit() {
    const w = root.clientWidth - (rail ? rail.offsetWidth : 0), h = window.innerHeight;
    const s = Math.min(w / W, h / H);
    if (Math.abs(s - scale) < 1e-4) return;
    scale = s;
    root.style.setProperty("--s", s.toFixed(5));
    document.dispatchEvent(new CustomEvent("slidescale", { detail: { scale: s } }));
  }
  fit();
  addEventListener("resize", fit);

  // ---- tabs: one panel shown at a time ----
  for (const box of document.querySelectorAll("[data-tabs]")) {
    const tabs = [...box.querySelectorAll('[role="tab"]')];
    // Selecting a tab also tells the page's own code ("tabchange"), so a drawing shared by the tabs can follow.
    const select = (tab) => {
      tabs.forEach((t) => {
        t.setAttribute("aria-selected", t === tab ? "true" : "false");
        document.getElementById(t.getAttribute("aria-controls")).hidden = t !== tab;
      });
      box.dispatchEvent(new CustomEvent("tabchange", { bubbles: true, detail: { tab: tab.id } }));
    };
    tabs.forEach((t) => t.addEventListener("click", () => select(t)));
    select(tabs.find((t) => t.getAttribute("aria-selected") === "true") || tabs[0]);
    // A "Next" button at the end of a panel opens the following tab.
    box.addEventListener("click", (e) => {
      const next = e.target.closest("[data-next]");
      if (next) select(document.getElementById(next.dataset.next));
    });
  }

  // ---- current page: the one crossing the middle of the window ----
  const frames = [...document.querySelectorAll(".frame[id]")];
  const links = [...document.querySelectorAll(".rail a[href^='#']")];
  let current = -1;
  const io = new IntersectionObserver((entries) => {
    for (const e of entries) {
      if (!e.isIntersecting) continue;
      const i = frames.indexOf(e.target);
      if (i === current) continue;
      current = i;
      links.forEach((a) => (a.hash === `#${e.target.id}` ? a.setAttribute("aria-current", "true") : a.removeAttribute("aria-current")));
      document.dispatchEvent(new CustomEvent("framechange", { detail: { id: e.target.id } }));
    }
  }, { rootMargin: "-50% 0px -50% 0px" });
  frames.forEach((f) => io.observe(f));
})();
