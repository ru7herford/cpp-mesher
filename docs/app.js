// The cpp-mesher page. Everything drawn here comes from docs/data/, written by scripts/export_page.py from the
// C++ program's own output: steps.json (every change to the mesh), report.json, progress.json (quality after
// every refinement step), checks.json (the three test programs), limits.json and triangle.json.
// The browser only replays and draws; it never meshes.
// The page is a deck of 1600 x 900 slides scaled by deck.js (CSS transform). Pointer positions are turned into
// slide pixels with localPoint(), and canvases are redrawn when the scale changes so they stay sharp.

const STAGES = [
  { id: "input", label: "1. Points on edges" },
  { id: "insert", label: "2. Connect points" },
  { id: "recover", label: "3. Put edges back" },
  { id: "exterior", label: "4. Remove outside" },
  { id: "legalize", label: "4. Tidy up" },
  { id: "refine", label: "5. Fix thin triangles", heading: true },
  { id: "split", label: "Split an edge", sub: true },
  { id: "circumcentre", label: "Add a point", sub: true },
  { id: "blocked", label: "Too close: split", sub: true },
  { id: "done", label: "Done" },
];
const REASON_STAGE = {
  encroached: "split", small_angle: "circumcentre", too_long: "circumcentre",
  centre_encroaches: "blocked", centre_outside: "blocked", centre_on_segment: "blocked",
};
const REPLAY_SECONDS = 25; // a whole run plays in about this long, whatever its size
const SETTLE_MS = 500;     // a frame must stay current this long before its replay starts by itself

const $ = (id) => document.getElementById(id);
const css = (name) => getComputedStyle(document.documentElement).getPropertyValue(name).trim();
const TOKENS = ["--panel", "--ink", "--ink-3", "--model", "--baseline"];
const tone = Object.fromEntries(TOKENS.map((name) => [name, css(name)]));
const ICONS = {
  play: '<svg viewBox="0 0 14 14" aria-hidden="true"><path d="M3 1.5v11l9.5-5.5z"/></svg>',
  pause: '<svg viewBox="0 0 14 14" aria-hidden="true"><path d="M2.5 1.5h3.2v11H2.5zM8.3 1.5h3.2v11H8.3z"/></svg>',
  zoomIn: '<svg viewBox="0 0 20 20" aria-hidden="true"><path d="M4 10h12M10 4v12"/></svg>',
  zoomOut: '<svg viewBox="0 0 20 20" aria-hidden="true"><path d="M4 10h12"/></svg>',
};
const key = (t) => t.join(",");
const fmt = (x, digits = 1) => x.toFixed(digits);
const int = (x) => x.toLocaleString("en-GB");
const plural = (n, word) => `${int(n)} ${word}${n === 1 ? "" : "s"}`;

/** A pointer event's position inside `el`, in the element's own (unscaled slide) pixels. The slide is scaled
 *  by a CSS transform, so the on-screen rectangle is divided back by its ratio to the layout size. */
function localPoint(e, el) {
  const r = el.getBoundingClientRect();
  return [((e.clientX - r.left) * el.offsetWidth) / r.width, ((e.clientY - r.top) * el.offsetHeight) / r.height];
}

async function json(path) {
  const r = await fetch(path);
  if (!r.ok) throw new Error(`${path}: ${r.status}`);
  return r.json();
}

// ---- colour ---------------------------------------------------------------------

const rgb = (h) => [1, 3, 5].map((i) => parseInt(h.slice(i, i + 2), 16));
function mix(a, b, t) {
  const [p, q] = [rgb(a), rgb(b)];
  return `rgb(${p.map((v, i) => Math.round(v + (q[i] - v) * t)).join(",")})`;
}
const alpha = (h, a) => `rgba(${rgb(h).join(",")},${a})`;

/** Quality fill: the model blue mixed into the panel, pale just above the target and darkest at 60°, in the
 *  same 5° classes as the histograms (progress.json bin_deg), so bars and tiles share exact colours. */
function palette(target, binDeg) {
  const paper = tone["--panel"], model = tone["--model"];
  const classes = Math.max(1, Math.ceil((60 - target) / binDeg));
  const ramp = Array.from({ length: classes }, (_, k) => mix(paper, model, 0.08 + (0.34 * k) / Math.max(1, classes - 1)));
  return {
    ramp, paper,
    edge: alpha(tone["--ink"], 0.3),
    bad: mix(paper, tone["--baseline"], 0.45),
    fill: (angle) => (angle < target ? null : ramp[Math.min(classes - 1, Math.floor((angle - target) / binDeg))]),
    ink: tone["--ink"], ink3: tone["--ink-3"], scaffold: alpha(tone["--ink"], 0.3), model, baseline: tone["--baseline"],
  };
}

// ---- geometry ---------------------------------------------------------------

/** The three angles in degrees, with the formula of src/quality.cpp:triangle_angles. */
function angles(p, q, r) {
  return [[p, q, r], [q, r, p], [r, p, q]].map(([c, n, m]) => {
    const ux = n[0] - c[0], uy = n[1] - c[1], vx = m[0] - c[0], vy = m[1] - c[1];
    return (Math.atan2(Math.abs(ux * vy - uy * vx), ux * vx + uy * vy) * 180) / Math.PI;
  });
}

function contains(a, b, c, p) {
  const cross = (u, v) => (v[0] - u[0]) * (p[1] - u[1]) - (v[1] - u[1]) * (p[0] - u[0]);
  const d1 = cross(a, b), d2 = cross(b, c), d3 = cross(c, a);
  return !((d1 < 0 || d2 < 0 || d3 < 0) && (d1 > 0 || d2 > 0 || d3 > 0));
}

function circumcircle(p, q, r) {
  const d = 2 * (p[0] * (q[1] - r[1]) + q[0] * (r[1] - p[1]) + r[0] * (p[1] - q[1]));
  const s = (v) => v[0] * v[0] + v[1] * v[1];
  const x = (s(p) * (q[1] - r[1]) + s(q) * (r[1] - p[1]) + s(r) * (p[1] - q[1])) / d;
  const y = (s(p) * (r[0] - q[0]) + s(q) * (p[0] - r[0]) + s(r) * (q[0] - p[0])) / d;
  return { x, y, r: Math.hypot(p[0] - x, p[1] - y) };
}

// ---- timeline: steps.json -> one event per thing that happened ----------------

function buildTimeline(run, report, example) {
  const n = run.points.length;
  const buildPoints = [...run.points, ...run.build.enclosing];
  const finalPoints = [...run.points, ...run.steps.map((s) => s.point)];
  const { input } = run;
  const events = [];
  const add = (stage, fields) => events.push({ stage, removed: [], added: [], points: buildPoints, ...fields });

  add("input", { text: `<b>Step 1.</b> ${example.note} Its edges are marked by ${n} points.` });

  run.build.insertions.forEach((change, i) => add("insert", {
    ...change, point: i, circles: change.removed,
    text: `<b>Step 2.</b> Point ${i + 1} of ${n} goes in. The triangles around it are removed (dashed) and rebuilt so they all meet at the new point.`,
  }));
  if (run.build.recoveries.length === 0) {
    add("recover", { empty: true, text: "<b>Step 3.</b> Every edge of the shape is already in the mesh. Nothing to put back." });
  }
  for (const { segment, removed, added } of run.build.recoveries) {
    add("recover", { removed, added, segment,
      text: `<b>Step 3.</b> An edge of the shape (thick) was missing. Flipping ${plural(removed.length, "triangle")} that crossed it puts it back.` });
  }
  add("exterior", { ...run.build.exterior,
    text: `<b>Step 4.</b> ${plural(run.build.exterior.removed.length, "triangle")} lie outside the shape${input.holes.length ? " or in its hole" : ""}. They are deleted.` });
  const flips = run.build.legalize.removed.length;
  add("legalize", { ...run.build.legalize, empty: flips === 0,
    text: flips === 0 ? "<b>Step 4.</b> This is the first full mesh. Orange triangles are thin."
      : `<b>Step 4.</b> A few flips tidy up the triangles (${plural(flips, "triangle")} changed). This is the first full mesh. Orange triangles are thin.` });
  const built = events.length - 1; // index of the first full mesh (refinement step 0)

  for (const step of run.steps) {
    add(REASON_STAGE[step.reason], { removed: step.removed, added: step.added, points: finalPoints, point: step.index, step,
      text: explainStep(step, finalPoints, input) });
  }
  const r = report.refinement, corner = r.of_which_at_sharp_input_corners;
  add("done", { points: finalPoints,
    text: `<b>Done.</b> ${int(report.triangles)} triangles. The smallest angle is ${fmt(report.min_angle_deg, 1)}° (target ${input.min_angle_deg}°).
      ${corner ? `The ${plural(corner, "triangle")} under the target sit in the ${fmt(Math.min(...cornerAngles(input)), 0)}° corner of the shape itself. No point can fix those.` : ""}` });
  const exterior = events.findIndex((e) => e.stage === "exterior");
  return { events, n, finalPoints, exterior, built, steps: run.steps.length };
}

/** Interior angles of the input corners, in degrees. */
function cornerAngles(input) {
  const out = [];
  [input.outer, ...input.holes].forEach((loop, l) => {
    let area = 0;
    loop.forEach((p, i) => { const q = loop[(i + 1) % loop.length]; area += p[0] * q[1] - q[0] * p[1]; });
    const ring = (area > 0) === (l === 0) ? loop : [...loop].reverse();
    ring.forEach((c, k) => {
      const nx = ring[(k + 1) % ring.length], pv = ring[(k + ring.length - 1) % ring.length];
      const a = (Math.atan2((nx[0] - c[0]) * (pv[1] - c[1]) - (nx[1] - c[1]) * (pv[0] - c[0]), (nx[0] - c[0]) * (pv[0] - c[0]) + (nx[1] - c[1]) * (pv[1] - c[1])) * 180) / Math.PI;
      out.push(a < 0 ? a + 360 : a);
    });
  });
  return out;
}

function explainStep(step, points, input) {
  const tri = step.triangle && step.triangle.map((i) => points[i]);
  const where = step.shell ? "near the corner" : "in half";
  const smallest = () => Math.min(...angles(...tri));
  const why = () => (smallest() < input.min_angle_deg ? `has a ${fmt(smallest())}° angle` : "is too big");
  switch (step.reason) {
    case "encroached":
      return `<b>Step 5.</b> A point (ringed) is too close to an edge of the shape. That edge is split ${where} first.`;
    case "small_angle":
      return `<b>Step 5.</b> This triangle is thin: a ${fmt(smallest())}° angle. A new point goes in at its circumcentre (×), which removes it.`;
    case "too_long":
      return `<b>Step 5.</b> This triangle is too big for the size limit. A new point goes in at its circumcentre (×).`;
    case "centre_encroaches":
      return `<b>Step 5.</b> This triangle ${why()}, but its new point (×) would be too close to an edge of the shape. The edge is split ${where} instead.`;
    case "centre_outside":
      return `<b>Step 5.</b> This triangle ${why()}, but its new point (×) would be outside the shape. The edge in the way is split instead.`;
    case "centre_on_segment":
      return "<b>Step 5.</b> The new point would be exactly on an edge of the shape, so that edge is split there.";
    default:
      return "";
  }
}

// ---- replay state ----------------------------------------------------------------

class Replay {
  constructor(timeline) { this.timeline = timeline; this.reset(); }
  reset() {
    const n = this.timeline.n;
    this.live = new Map([[key([n, n + 1, n + 2]), [n, n + 1, n + 2]]]); // the enclosing triangle
    this.index = -1;
    this.pointCount = 0;
    this.cache = new Map();
  }
  /** Applies events up to and including `target`. Going back replays from the start (cheap). */
  seek(target) {
    if (target < this.index) this.reset();
    while (this.index < target) {
      const event = this.timeline.events[++this.index];
      for (const t of event.removed) this.live.delete(key(t));
      for (const t of event.added) this.live.set(key(t), t);
      if (event.point !== undefined) this.pointCount = event.point + 1;
    }
  }
  get event() { return this.timeline.events[this.index]; }
  /** True once the outside triangles are gone, so angles mean something. */
  get clipped() { return this.index >= this.timeline.exterior; }
  /** Smallest angle of a live triangle. Cached by vertex ids only once the outside is gone:
   *  before that, ids n..n+2 are the enclosing triangle's corners rather than refinement points. */
  smallest(t) {
    const pts = this.event.points;
    if (!this.clipped) return Math.min(...angles(pts[t[0]], pts[t[1]], pts[t[2]]));
    const k = key(t);
    let a = this.cache.get(k);
    if (a === undefined) this.cache.set(k, (a = Math.min(...angles(pts[t[0]], pts[t[1]], pts[t[2]]))));
    return a;
  }
  /** Refinement step shown (0 = first full mesh), or -1 while the mesh is still being built. */
  get step() { return this.index < this.timeline.built ? -1 : Math.min(this.index - this.timeline.built, this.timeline.steps); }
}

// ---- drawing ---------------------------------------------------------------------

function makeView(canvas, box) {
  const rect = canvas.getBoundingClientRect();
  const dpr = Math.min(window.devicePixelRatio || 1, 2);
  canvas.width = Math.max(1, Math.round(rect.width * dpr));
  canvas.height = Math.max(1, Math.round(rect.height * dpr));
  const [x0, y0, x1, y1] = box;
  const pad = (10 * canvas.width) / Math.max(1, canvas.offsetWidth); // 10 stage pixels, so the drawing is the same at every scale
  const scale = Math.min((canvas.width - 2 * pad) / (x1 - x0), (canvas.height - 2 * pad) / (y1 - y0));
  const ox = (canvas.width - scale * (x1 - x0)) / 2, oy = (canvas.height - scale * (y1 - y0)) / 2;
  return {
    dpr, scale,
    x: (p) => ox + (p[0] - x0) * scale,
    y: (p) => canvas.height - oy - (p[1] - y0) * scale,
    unproject: (px, py) => [x0 + (px - ox) / scale, y0 + (canvas.height - oy - py) / scale],
  };
}

/** Draws the replay's current mesh. `plain` leaves out what the current event changed and why. Two extras for
 *  the "before" pictures of page 4: `segment`, an edge of the shape to mark, and `shade`, triangles to grey out. */
function draw(canvas, view, replay, { palette: colour, loops, n, probe = null, plain = false, segment = null, shade = [] }) {
  const ctx = canvas.getContext("2d");
  const { dpr } = view;
  const event = replay.event;
  const pts = event.points;
  const path = (t) => {
    ctx.beginPath();
    ctx.moveTo(view.x(pts[t[0]]), view.y(pts[t[0]]));
    ctx.lineTo(view.x(pts[t[1]]), view.y(pts[t[1]]));
    ctx.lineTo(view.x(pts[t[2]]), view.y(pts[t[2]]));
    ctx.closePath();
  };
  const circle = (c) => {
    ctx.beginPath();
    ctx.setLineDash([4 * dpr, 4 * dpr]);
    ctx.arc(view.x([c.x, c.y]), view.y([c.x, c.y]), c.r * view.scale, 0, 2 * Math.PI);
    ctx.stroke();
    ctx.setLineDash([]);
  };
  const dot = (p, radius, fill) => {
    ctx.beginPath();
    ctx.arc(view.x(p), view.y(p), radius * dpr, 0, 2 * Math.PI);
    ctx.fillStyle = fill;
    ctx.fill();
  };
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.lineJoin = "round";
  const filled = replay.clipped;

  ctx.fillStyle = alpha(colour.ink, 0.12);
  for (const t of shade) { path(t); ctx.fill(); }
  ctx.lineWidth = (filled ? 0.6 : 0.75) * dpr;
  for (const t of replay.live.values()) {
    path(t);
    if (filled) {
      ctx.fillStyle = colour.fill(replay.smallest(t)) ?? colour.bad;
      ctx.fill();
      ctx.strokeStyle = colour.edge;
    } else ctx.strokeStyle = colour.scaffold;
    ctx.stroke();
  }

  if (!plain) {
    // What this event changed: new triangles outlined (tinted while there is no fill), removed ones dashed.
    if (!filled) {
      ctx.globalAlpha = 0.2;
      ctx.fillStyle = colour.model;
      for (const t of event.added) { path(t); ctx.fill(); }
      ctx.globalAlpha = 1;
    }
    ctx.strokeStyle = filled ? colour.ink : colour.model;
    ctx.lineWidth = 1.5 * dpr;
    for (const t of event.added) { path(t); ctx.stroke(); }
    ctx.strokeStyle = colour.ink3;
    ctx.setLineDash([3 * dpr, 3 * dpr]);
    ctx.lineWidth = dpr;
    const MAX_DASHED = 400; // past this many (the big exterior sweep) the dashes would bury the drawing
    if (event.stage !== "exterior" || event.removed.length < MAX_DASHED) for (const t of event.removed) { path(t); ctx.stroke(); }
    ctx.setLineDash([]);
  }

  // The input boundary, always on top.
  ctx.strokeStyle = colour.ink;
  ctx.lineWidth = 2 * dpr;
  for (const loop of loops) {
    ctx.beginPath();
    loop.forEach((p, i) => (i ? ctx.lineTo(view.x(p), view.y(p)) : ctx.moveTo(view.x(p), view.y(p))));
    ctx.closePath();
    ctx.stroke();
  }

  // The cause of the step, in the colour of its kind: circumcentres in model blue, splits in baseline orange.
  const step = plain ? null : event.step;
  const edge = plain ? segment : step?.segment || event.segment;
  if (edge) {
    const [a, b] = edge.map((i) => pts[i]);
    ctx.strokeStyle = colour.baseline;
    ctx.lineWidth = 4 * dpr;
    ctx.beginPath();
    ctx.moveTo(view.x(a), view.y(a));
    ctx.lineTo(view.x(b), view.y(b));
    ctx.stroke();
  }
  if (!plain) {
    ctx.strokeStyle = step && REASON_STAGE[step.reason] === "circumcentre" ? colour.model : colour.baseline;
    ctx.lineWidth = 1.5 * dpr;
    if (event.circles) { // step 2: the circles of the triangles a new point removes, in grey (nothing is wrong)
      ctx.strokeStyle = colour.ink3;
      for (const t of event.circles) circle(circumcircle(pts[t[0]], pts[t[1]], pts[t[2]]));
    }
    if (step?.triangle) {
      path(step.triangle);
      ctx.globalAlpha = 0.35;
      ctx.fillStyle = ctx.strokeStyle;
      ctx.fill();
      ctx.globalAlpha = 1;
      ctx.stroke();
      circle(circumcircle(...step.triangle.map((i) => pts[i])));
    }
    if (step?.segment) {
      const [a, b] = step.segment.map((i) => pts[i]);
      ctx.strokeStyle = colour.baseline;
      circle({ x: (a[0] + b[0]) / 2, y: (a[1] + b[1]) / 2, r: Math.hypot(a[0] - b[0], a[1] - b[1]) / 2 });
    }
    if (step?.encroacher !== undefined) {
      const p = pts[step.encroacher];
      ctx.beginPath();
      ctx.arc(view.x(p), view.y(p), 6 * dpr, 0, 2 * Math.PI);
      ctx.stroke();
    }
    if (step?.centre) {
      const s = 5 * dpr, x = view.x(step.centre), y = view.y(step.centre);
      ctx.beginPath();
      ctx.moveTo(x - s, y - s); ctx.lineTo(x + s, y + s);
      ctx.moveTo(x - s, y + s); ctx.lineTo(x + s, y - s);
      ctx.stroke();
    }
  }

  if (probe) {
    path(probe);
    ctx.strokeStyle = colour.ink;
    ctx.lineWidth = 2.5 * dpr;
    ctx.stroke();
  }

  if (!filled) {
    const shown = event.stage === "input" ? n : replay.pointCount;
    for (let i = 0; i < Math.min(shown, pts.length); i++) dot(pts[i], i < n ? 1.6 : 1.1, colour.ink);
  }
  if (!plain && event.point !== undefined) {
    dot(pts[event.point], 4.5, colour.paper);
    dot(pts[event.point], 3.5, colour.ink);
  }
}

// ---- charts (SVG) ------------------------------------------------------------------

const NS = "http://www.w3.org/2000/svg";
function svgOf(host, height) {
  const width = Math.max(220, host.clientWidth || 300);
  host.innerHTML = `<svg viewBox="0 0 ${width} ${height}" width="${width}" height="${height}"></svg>`;
  return [host.firstChild, width];
}
function el(svg, name, attributes, text) {
  const node = document.createElementNS(NS, name);
  for (const [k, v] of Object.entries(attributes)) node.setAttribute(k, v);
  if (text !== undefined) node.textContent = text;
  svg.appendChild(node);
  return node;
}

/** Thin triangles left after each new point of step 5, with a cursor at the current point. */
function drawQuality(host, progress, step) {
  const H = 240, [svg, W] = svgOf(host, H);
  const m = { l: 40, r: 12, t: 16, b: 30 };
  const N = progress.below.length - 1, max = Math.max(1, ...progress.below);
  const X = (s) => m.l + (s / Math.max(1, N)) * (W - m.l - m.r);
  const Y = (v) => m.t + (1 - v / max) * (H - m.t - m.b);
  el(svg, "line", { x1: m.l, x2: W - m.r, y1: Y(0), y2: Y(0), class: "axis" });
  for (const v of [0, max]) el(svg, "text", { x: m.l - 8, y: Y(v) + 6, "text-anchor": "end" }, v);
  el(svg, "polyline", { points: progress.below.map((v, s) => `${fmt(X(s), 1)},${fmt(Y(v), 1)}`).join(" "), class: "series" });
  if (step >= 0) {
    const v = progress.below[step];
    el(svg, "line", { x1: X(step), x2: X(step), y1: m.t, y2: Y(0), class: "cursor" });
    el(svg, "circle", { cx: X(step), cy: Y(v), r: 5, class: "handle" });
    const right = X(step) > W - 60;
    el(svg, "text", { x: X(step) + (right ? -10 : 10), y: Math.max(m.t + 14, Y(v) - 10), "text-anchor": right ? "end" : "start", class: "value" }, v);
  }
  for (const s of [0, N]) el(svg, "text", { x: X(s), y: H - 6, "text-anchor": s === 0 ? "start" : "end" }, s === 0 ? "first mesh" : `${int(s)} points`);

  // Where a pointer event lands on the step axis (the chart is redrawn every step, so this is handed back).
  // localPoint gives stage pixels, and the SVG is drawn at its layout width W, so they match directly.
  return (e) => {
    const [x] = localPoint(e, host);
    return Math.round(Math.max(0, Math.min(N, ((x - m.l) / (W - m.l - m.r)) * N)));
  };
}

// ---- page ------------------------------------------------------------------------

const state = { listing: null, runs: new Map(), id: null };

async function load(id) {
  if (!state.runs.has(id)) {
    const [run, report, progress] = await Promise.all(["steps", "report", "progress"].map((f) => json(`data/${id}/${f}.json`)));
    const example = state.listing.examples.find((e) => e.id === id);
    const timeline = buildTimeline(run, report, example);
    state.runs.set(id, { run, report, progress, example, timeline, loops: [run.input.outer, ...run.input.holes],
      box: bounds(run.input.outer), palette: palette(run.input.min_angle_deg, progress.bin_deg) });
  }
  return state.runs.get(id);
}

function bounds(points) {
  const xs = points.map((p) => p[0]), ys = points.map((p) => p[1]);
  return [Math.min(...xs), Math.min(...ys), Math.max(...xs), Math.max(...ys)];
}

// Page 2: a close-up of the airfoil's first full mesh, before any thin triangle is fixed, with the thinnest
// triangle outlined and labelled with its angle.
function drawFirst(cur) {
  const replay = new Replay(cur.timeline);
  replay.seek(cur.timeline.built);
  const canvas = $("first-mesh");
  const view = makeView(canvas, cur.example.focus || cur.box);
  let thinnest = null, angle = Infinity;
  for (const t of replay.live.values()) if (replay.smallest(t) < angle) { angle = replay.smallest(t); thinnest = t; }
  draw(canvas, view, replay, { palette: cur.palette, loops: cur.loops, n: cur.timeline.n, probe: thinnest, plain: true });

  // The label sits just above the middle of the triangle's shortest edge (its visible short end), kept inside the picture.
  const pts = thinnest.map((i) => replay.event.points[i]);
  const px = (s) => (s * canvas.width) / canvas.offsetWidth; // slide pixels -> canvas pixels
  const clamp = (v, lo, hi) => Math.max(lo, Math.min(hi, v));
  const [a, b] = [[0, 1], [1, 2], [2, 0]].map(([i, j]) => [pts[i], pts[j]])
    .sort((e, f) => Math.hypot(e[0][0] - e[1][0], e[0][1] - e[1][1]) - Math.hypot(f[0][0] - f[1][0], f[0][1] - f[1][1]))[0];
  const x = clamp((view.x(a) + view.x(b)) / 2 + px(10), px(16), canvas.width - px(240));
  const y = clamp((view.y(a) + view.y(b)) / 2 - px(24), px(30), canvas.height - px(16));
  const ctx = canvas.getContext("2d");
  ctx.font = `600 ${px(20)}px ${css("--font")}`;
  const text = `thinnest: ${fmt(angle, 3)}°`;
  ctx.fillStyle = tone["--panel"];
  ctx.fillRect(x - px(6), y - px(22), ctx.measureText(text).width + px(12), px(30));
  ctx.fillStyle = tone["--ink"];
  ctx.fillText(text, x, y);
}

// Page 4: the five steps on the slot, one tab each. Each tab has its own list of moments to step through:
// its points, its events, or a "before" picture followed by what the step did.
const how = { cur: null, replay: null, tab: "how-1", at: {}, moments: {} };
const stripStep = (text) => text.replace(/^<b>Step \d\.<\/b> /, "");

function howMoments(cur) {
  const { events, n, built } = cur.timeline;
  const at = (stage) => events.map((e, i) => (e.stage === stage && !e.empty ? i : -1)).filter((i) => i >= 0);
  const inserts = at("insert"), recovers = at("recover"), exterior = at("exterior")[0];
  const missing = events[recovers[0]];
  const last = (i) => ({ index: i, text: stripStep(events[i].text) });
  const steps = cur.run.steps.length, p = cur.progress;
  return {
    "how-1": Array.from({ length: n }, (_, k) => ({ points: k + 1, label: `point ${k + 1} of ${n}`,
      text: k + 1 < n ? `Point ${k + 1}: a corner of the shape.` : `All ${n} points. The straight lines between them are the shape's edges.` })),
    "how-2": inserts.map((i, k) => ({ ...last(i), label: `point ${k + 1} of ${n}` })),
    "how-3": [{ index: recovers[0] - 1, plain: true, segment: missing.segment, label: "before",
      text: `The orange edge belongs to the shape, but ${plural(missing.removed.length, "triangle")} cross it.` },
    ...recovers.map((i) => ({ ...last(i), label: "after" }))],
    "how-4": [{ index: exterior - 1, plain: true, shade: events[exterior].removed, label: "before",
      text: `The ${events[exterior].removed.length} grey triangles lie outside the shape.` },
    { ...last(exterior), label: "outside removed" }, { ...last(built), label: "tidied" }],
    "how-5": [{ index: built, plain: true, label: "first mesh", text: `The first full mesh: ${p.below[0]} of ${p.triangles[0]} triangles are thin (orange).` },
      ...cur.run.steps.map((_, k) => ({ ...last(built + 1 + k), label: `fix ${k + 1} of ${steps}` }))],
  };
}

/** Step 1's picture: the whole shape faintly, and its first `count` points joined by its edges. */
function drawPoints(canvas, view, cur, count) {
  const ctx = canvas.getContext("2d"), { dpr } = view;
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.strokeStyle = alpha(tone["--ink"], 0.2);
  ctx.lineWidth = 2 * dpr;
  for (const loop of cur.loops) {
    ctx.beginPath();
    loop.forEach((p, i) => (i ? ctx.lineTo(view.x(p), view.y(p)) : ctx.moveTo(view.x(p), view.y(p))));
    ctx.closePath();
    ctx.stroke();
  }
  let shown = 0;
  for (const loop of cur.loops) {
    const k = Math.min(loop.length, count - shown);
    shown += loop.length;
    if (k <= 0) continue;
    ctx.strokeStyle = tone["--ink"];
    ctx.lineWidth = 2 * dpr;
    ctx.beginPath();
    loop.slice(0, k).forEach((p, i) => (i ? ctx.lineTo(view.x(p), view.y(p)) : ctx.moveTo(view.x(p), view.y(p))));
    if (k === loop.length) ctx.closePath();
    ctx.stroke();
    loop.slice(0, k).forEach((p, i) => {
      ctx.beginPath();
      ctx.arc(view.x(p), view.y(p), (i === k - 1 && k < loop.length ? 6 : 4) * dpr, 0, 2 * Math.PI);
      ctx.fillStyle = tone["--ink"];
      ctx.fill();
    });
  }
}

function drawHow() {
  const cur = how.cur, list = how.moments[how.tab], k = how.at[how.tab], m = list[k];
  const canvas = $("how-mesh"), view = makeView(canvas, cur.box);
  if (m.points) drawPoints(canvas, view, cur, m.points);
  else {
    how.replay.seek(m.index);
    draw(canvas, view, how.replay, { palette: cur.palette, loops: cur.loops, n: cur.timeline.n, plain: m.plain, segment: m.segment, shade: m.shade });
  }
  $("how-pos").textContent = m.label;
  $("how-note").innerHTML = m.text;
  $("how-prev").disabled = k === 0;
  $("how-next").disabled = k === list.length - 1;
}

function initHow(cur) {
  how.cur = cur;
  how.replay = new Replay(cur.timeline);
  how.moments = howMoments(cur);
  for (const tab of Object.keys(how.moments)) how.at[tab] = 0;
  how.tab = document.querySelector('#how [role="tab"][aria-selected="true"]').id;
  drawHow();
}

// Page 5: the replay of the chosen shape, with its results.
const watch = { replay: null, playing: false, timer: 0, away: false, rate: 1, stepAt: null, view: null, zoomed: false, probe: null };

function renderStages() {
  const { events } = state.runs.get(state.id).timeline;
  $("stages").innerHTML = STAGES.map((s) => {
    if (s.heading) return `<li class="group">${s.label}</li>`;
    const has = s.id === "input" || s.id === "done" || events.some((e) => e.stage === s.id && !e.empty);
    return `<li class="${s.sub ? "sub" : ""}"><button type="button" data-stage="${s.id}" class="${has ? "" : "empty"}">
      <span>${s.label}</span><span class="meta"></span></button></li>`;
  }).join("");
}

/** Work done by each stage in events 0..upto: points, edges, or (for "Remove outside") triangles removed. */
function stageCounts(events, upto) {
  const counts = {};
  for (let i = 0; i <= upto; i++) {
    const e = events[i];
    if (!e.empty) counts[e.stage] = (counts[e.stage] || 0) + (e.stage === "exterior" ? e.removed.length : 1);
  }
  return counts;
}

function renderWatch() {
  const cur = state.runs.get(state.id);
  const { replay } = watch;
  const { timeline, progress } = cur;
  const event = replay.event;
  const canvas = $("watch-mesh");
  if (watch.probe && !(replay.clipped && replay.live.has(key(watch.probe)))) watch.probe = null;
  watch.view = makeView(canvas, watch.zoomed ? cur.example.focus : cur.box);
  draw(canvas, watch.view, replay, { palette: cur.palette, loops: cur.loops, n: timeline.n, probe: watch.probe });
  renderProbe();

  $("explain").innerHTML = event.text;
  $("explain").setAttribute("aria-live", watch.playing ? "off" : "polite");
  const last = timeline.events.length - 1, s = replay.step;
  $("scrub").value = replay.index;
  $("readout").textContent = s < 0 ? "steps 1–4" : `${int(s)} / ${int(timeline.steps)}`;
  $("prev").disabled = replay.index === 0;
  $("next").disabled = replay.index === last;

  const counts = stageCounts(timeline.events, replay.index), totals = stageCounts(timeline.events, last);
  for (const b of $("stages").querySelectorAll("button")) {
    const id = b.dataset.stage;
    if (id === event.stage) b.setAttribute("aria-current", "step"); else b.removeAttribute("aria-current");
    b.querySelector(".meta").textContent = id === "input" || id === "done" || !totals[id] ? "" : `${int(counts[id] ?? 0)} / ${int(totals[id])}`;
  }

  watch.stepAt = drawQuality($("quality"), progress, s);
  setPlay();
}

/** Under the mesh: the angles of the triangle under the pointer, or how to get them. Shown once the outside is
 *  gone, since before that the triangles are only scaffolding. */
function renderProbe() {
  const probe = $("probe"), clipped = watch.replay.clipped;
  probe.hidden = !clipped;
  $("watch-mesh").classList.toggle("probing", clipped);
  probe.classList.toggle("hint", !watch.probe);
  if (!watch.probe) {
    probe.textContent = matchMedia("(hover: hover)").matches ? "Point at a triangle to read its angles." : "Tap a triangle to read its angles.";
    return;
  }
  const target = state.runs.get(state.id).run.input.min_angle_deg;
  const angs = angles(...watch.probe.map((i) => watch.replay.event.points[i])).sort((a, b) => a - b);
  probe.innerHTML = `This triangle: <b>${fmt(angs[0])}°</b>, ${fmt(angs[1])}°, ${fmt(angs[2])}°. ${angs[0] < target ? "Thin." : "Not thin."}`;
}

function pickTriangle(e) {
  if (!watch.replay.clipped) return null;
  const canvas = $("watch-mesh");
  const [x, y] = localPoint(e, canvas); // stage pixels; the drawing buffer has its own size
  const m = watch.view.unproject((x * canvas.width) / canvas.offsetWidth, (y * canvas.height) / canvas.offsetHeight);
  const pts = watch.replay.event.points;
  for (const t of watch.replay.live.values()) if (contains(pts[t[0]], pts[t[1]], pts[t[2]], m)) return t;
  return null;
}

/** The zoom button shows its next action: + to zoom in, − to zoom back out. */
function setZoom() {
  $("zoom").innerHTML = watch.zoomed ? ICONS.zoomOut : ICONS.zoomIn;
  $("zoom").setAttribute("aria-label", watch.zoomed ? "Zoom out" : "Zoom in");
}

/** The results box for the chosen shape, and, for the sharp corner, why its smallest angle is under the target. */
function renderResults(cur) {
  const { example, report } = cur;
  $("results").innerHTML = `
    <dt>Triangles</dt><dd>${int(example.triangles)}</dd>
    <dt>Smallest angle</dt><dd>${fmt(example.min_angle_deg, 1)}°</dd>
    <dt>Points added</dt><dd>${int(example.steps)}</dd>
    <dt>Time</dt><dd>${fmt(example.time_ms, 1)} ms</dd>`;
  const corner = report.refinement.of_which_at_sharp_input_corners;
  $("results-note").textContent = corner
    ? `The ${plural(corner, "triangle")} under ${report.targets.min_angle_deg}° sit in the shape's own ${fmt(Math.min(...cornerAngles(cur.run.input)), 0)}° corner. No mesher can do better there.`
    : "";
}

function setPlay() {
  const cur = state.runs.get(state.id);
  const atEnd = watch.replay.index >= cur.timeline.events.length - 1;
  $("play").innerHTML = watch.playing ? `${ICONS.pause} Pause` : `${ICONS.play} ${atEnd ? "Replay" : "Play"}`;
}

function goTo(index) {
  const cur = state.runs.get(state.id);
  watch.replay.seek(Math.max(0, Math.min(cur.timeline.events.length - 1, index)));
  renderWatch();
}

function play() {
  const cur = state.runs.get(state.id);
  // From the end, replay the refinement from the first full mesh; the construction stages are a click away.
  if (watch.replay.index >= cur.timeline.events.length - 1) goTo(cur.timeline.built);
  watch.playing = true;
  let last = performance.now(), carry = 0;
  const tick = (now) => {
    if (!watch.playing) return;
    carry += ((now - last) / 1000) * watch.rate;
    last = now;
    if (carry >= 1) {
      goTo(watch.replay.index + Math.floor(carry));
      carry -= Math.floor(carry);
    }
    if (watch.replay.index >= cur.timeline.events.length - 1) return stop();
    watch.timer = requestAnimationFrame(tick);
  };
  watch.timer = requestAnimationFrame(tick);
  setPlay();
}

function stop() {
  watch.playing = false;
  watch.away = false;
  cancelAnimationFrame(watch.timer);
  if (watch.replay) setPlay();
}

/** Shows a shape on page 5, at its finished mesh. */
async function choose(id) {
  stop();
  const cur = await load(id);
  state.id = id;
  for (const input of document.querySelectorAll('input[name="shape"]')) input.checked = input.value === id;
  watch.replay = new Replay(cur.timeline);
  watch.rate = Math.max(4, cur.timeline.events.length / REPLAY_SECONDS);
  watch.zoomed = false;
  watch.probe = null;
  $("zoom").hidden = !cur.example.focus;
  setZoom();
  $("scrub").max = cur.timeline.events.length - 1;
  renderStages();
  renderResults(cur);
  goTo(cur.timeline.events.length - 1);
}

/** Numbers written into the text, read from the data. */
function renderFacts(airfoil, checks, limits, triangle) {
  const p = airfoil.progress;
  const targets = [...new Set(limits.runs.map((r) => r.target_deg))];
  const best = Math.max(...targets.filter((t) => limits.runs.filter((r) => r.target_deg === t).every((r) => r.met)));
  const ratios = triangle.examples.map((r) => r.triangle_ratio);
  const facts = {
    "airfoil.start_triangles": p.triangles[0], "airfoil.start_below": p.below[0], "airfoil.start_min": fmt(p.min_angle[0], 3),
    "limits.best": best, "limits.max_points": int(limits.max_points),
    "checks.naive_wrong": int(checks.predicates.orient_grid_naive_wrong), "checks.grid": int(checks.predicates.orient_grid),
    "checks.cases": int(checks.predicates.cases), "checks.fuzz": int(checks.fuzz.meshes),
    "triangle.range": `${fmt(Math.min(...ratios), 1)} to ${fmt(Math.max(...ratios), 1)}`,
  };
  for (const node of document.querySelectorAll("[data-fact]")) node.textContent = facts[node.dataset.fact];
}

function wire() {
  // Page 4: the stepper walks through the open tab's moments; each tab remembers where it was.
  document.addEventListener("tabchange", (e) => { if (how.cur && e.detail.tab in how.moments) { how.tab = e.detail.tab; drawHow(); } });
  const walk = (d) => { how.at[how.tab] = Math.max(0, Math.min(how.moments[how.tab].length - 1, how.at[how.tab] + d)); drawHow(); };
  $("how-prev").addEventListener("click", () => walk(-1));
  $("how-next").addEventListener("click", () => walk(1));

  // Page 5.
  $("shapes").addEventListener("change", (e) => choose(e.target.value));
  $("zoom").addEventListener("click", () => { watch.zoomed = !watch.zoomed; setZoom(); renderWatch(); });
  const mesh = $("watch-mesh");
  const setProbe = (t) => {
    if ((t && key(t)) === (watch.probe && key(watch.probe))) return;
    watch.probe = t;
    renderWatch();
  };
  mesh.addEventListener("pointermove", (e) => { if (e.pointerType === "mouse") setProbe(pickTriangle(e)); });
  mesh.addEventListener("pointerleave", () => setProbe(null));
  mesh.addEventListener("click", (e) => setProbe(pickTriangle(e)));

  const step = (d) => { stop(); goTo(watch.replay.index + d); };
  $("play").addEventListener("click", () => (watch.playing ? stop() : play()));
  $("prev").addEventListener("click", () => step(-1));
  $("next").addEventListener("click", () => step(1));
  $("scrub").addEventListener("input", (e) => { const to = Number(e.target.value); stop(); goTo(to); });
  $("stages").addEventListener("click", (e) => {
    const button = e.target.closest("button");
    if (!button) return;
    stop();
    // The next event of this stage after the current one, wrapping around: repeated clicks walk through them.
    const { events } = state.runs.get(state.id).timeline;
    const here = watch.replay.index;
    const matches = events.map((ev, i) => (ev.stage === button.dataset.stage ? i : -1)).filter((i) => i >= 0);
    if (matches.length) goTo(matches.find((i) => i > here) ?? matches[0]);
  });
  // Click or drag on the quality chart to jump to that refinement step.
  const chart = $("quality");
  const jump = (e) => { stop(); goTo(state.runs.get(state.id).timeline.built + watch.stepAt(e)); };
  chart.addEventListener("pointerdown", (e) => {
    if (e.button !== 0 || !watch.stepAt) return;
    chart.setPointerCapture(e.pointerId);
    jump(e);
  });
  chart.addEventListener("pointermove", (e) => { if (chart.hasPointerCapture(e.pointerId)) jump(e); });
  $("replay").addEventListener("keydown", (e) => {
    if (e.target.matches("input")) return; // the slider and radios keep their own arrow keys
    if (e.key === " " && !e.target.closest("button")) { e.preventDefault(); watch.playing ? stop() : play(); }
    if (e.key === "ArrowLeft" || e.key === "ArrowRight") { e.preventDefault(); step(e.key === "ArrowLeft" ? -1 : 1); }
  });

  // Arriving on page 5 plays the run: from the first mesh if it had finished, or on from where it was paused by
  // leaving. A pause you chose stays paused (stop() clears `away`). Leaving the page pauses it.
  // The cover's button always starts the airfoil from its first mesh.
  let inWatch = false, settle = 0;
  const pauseAway = () => { if (watch.playing) { stop(); watch.away = true; } };
  const arrive = () => {
    const atEnd = watch.replay.index >= state.runs.get(state.id).timeline.events.length - 1;
    if (atEnd || watch.away) { watch.away = false; play(); }
  };
  $("start-watch").addEventListener("click", async () => {
    if (state.id !== "airfoil") await choose("airfoil");
    stop();
    goTo(state.runs.get("airfoil").timeline.built);
    watch.away = true;
  });
  document.addEventListener("framechange", (e) => {
    inWatch = e.detail.id === "watch";
    clearTimeout(settle);
    if (!inWatch) return pauseAway();
    settle = setTimeout(arrive, SETTLE_MS);
  });
  document.addEventListener("visibilitychange", () => (document.hidden ? pauseAway() : inWatch && arrive()));

  // The layout inside a slide never changes; only its scale does. Redraw the canvases then, so their pixels
  // match the screen again (one redraw per animation frame at most).
  let pending = 0;
  document.addEventListener("slidescale", () => {
    cancelAnimationFrame(pending);
    pending = requestAnimationFrame(() => {
      renderWatch();
      drawHow();
      drawFirst(state.runs.get("airfoil"));
    });
  });
}

async function main() {
  const [listing, checks, limits, triangle] = await Promise.all(["examples", "checks", "limits", "triangle"].map((f) => json(`data/${f}.json`)));
  state.listing = listing;
  $("shapes").insertAdjacentHTML("beforeend", listing.examples.map((e) =>
    `<label><input type="radio" name="shape" value="${e.id}"><span>${e.title}</span></label>`).join(""));
  const [airfoil, slot] = await Promise.all([load("airfoil"), load("slot")]);
  renderFacts(airfoil, checks, limits, triangle);
  wire();
  drawFirst(airfoil);
  initHow(slot);
  await choose("airfoil");
}

main().catch((error) => {
  $("explain").textContent = `Could not load the recorded data (${error.message}).`;
  console.error(error);
});
