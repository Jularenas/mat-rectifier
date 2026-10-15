'use strict';
// Grid-based rectification: finds the mat's printed 1 cm grid intersections and fits
// image = distort(H · (u,v,1)), with u,v in cm, H a homography and a one-parameter
// division-model lens correction. Works with partial views and close-up (wide-angle) shots.
const Grid = (() => {
  const dist = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1]);

  // ---------- small linear algebra ----------
  function solve(A, b) {               // Gaussian elimination with partial pivoting (A: n×n array of rows)
    const n = b.length, M = A.map((r, i) => [...r, b[i]]);
    for (let c = 0; c < n; c++) {
      let p = c;
      for (let r = c + 1; r < n; r++) if (Math.abs(M[r][c]) > Math.abs(M[p][c])) p = r;
      if (Math.abs(M[p][c]) < 1e-12) return null;
      [M[c], M[p]] = [M[p], M[c]];
      for (let r = c + 1; r < n; r++) {
        const f = M[r][c] / M[c][c];
        if (f) for (let k = c; k <= n; k++) M[r][k] -= f * M[c][k];
      }
    }
    const x = new Array(n);
    for (let r = n - 1; r >= 0; r--) {
      let s = M[r][n];
      for (let k = r + 1; k < n; k++) s -= M[r][k] * x[k];
      x[r] = s / M[r][r];
    }
    return x;
  }
  const mul3 = (A, B) => A.map((_, i) => [0, 1, 2].map(j => A[i][0] * B[0][j] + A[i][1] * B[1][j] + A[i][2] * B[2][j]));
  function inv3(m) {
    const [a, b, c] = m[0], [d, e, f] = m[1], [g, h, i] = m[2];
    const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
    const det = a * A + b * B + c * C;
    return [[A / det, -(b * i - c * h) / det, (b * f - c * e) / det],
            [B / det, (a * i - c * g) / det, -(a * f - c * d) / det],
            [C / det, -(a * h - b * g) / det, (a * e - b * d) / det]];
  }
  const applyH = (H, x, y) => {
    const w = H[2][0] * x + H[2][1] * y + H[2][2];
    return [(H[0][0] * x + H[0][1] * y + H[0][2]) / w, (H[1][0] * x + H[1][1] * y + H[1][2]) / w];
  };
  function normT(pts) {                 // Hartley normalisation
    let mx = 0, my = 0;
    for (const p of pts) { mx += p[0]; my += p[1]; }
    mx /= pts.length; my /= pts.length;
    let md = 0;
    for (const p of pts) md += Math.hypot(p[0] - mx, p[1] - my);
    const s = Math.SQRT2 / (md / pts.length || 1);
    return [[s, 0, -s * mx], [0, s, -s * my], [0, 0, 1]];
  }
  // least-squares homography src→dst (normalised DLT, h33 = 1)
  function homography(src, dst) {
    if (src.length < 4) return null;
    const Ts = normT(src), Td = normT(dst);
    const ATA = Array.from({ length: 8 }, () => new Array(8).fill(0)), ATb = new Array(8).fill(0);
    const add = (row, rhs) => {
      for (let i = 0; i < 8; i++) { ATb[i] += row[i] * rhs; for (let j = 0; j < 8; j++) ATA[i][j] += row[i] * row[j]; }
    };
    for (let k = 0; k < src.length; k++) {
      const [u, v] = applyH(Ts, src[k][0], src[k][1]), [x, y] = applyH(Td, dst[k][0], dst[k][1]);
      add([u, v, 1, 0, 0, 0, -u * x, -v * x], x);
      add([0, 0, 0, u, v, 1, -u * y, -v * y], y);
    }
    const h = solve(ATA, ATb);
    if (!h) return null;
    const Hn = [[h[0], h[1], h[2]], [h[3], h[4], h[5]], [h[6], h[7], 1]];
    return mul3(mul3(inv3(Td), Hn), Ts);
  }

  // ---------- model: image = distort(H(u,v)) ----------
  // Division model: undistorted = c + (d − c) / (1 + λ·r_d²), r normalised by R (half-diagonal).
  function makeModel(H, lam, c, R, corr = null) {
    const Hi = inv3(H);
    const undistort = (x, y) => {
      const dx = (x - c[0]) / R, dy = (y - c[1]) / R, f = 1 / (1 + lam * (dx * dx + dy * dy));
      return [c[0] + dx * R * f, c[1] + dy * R * f];
    };
    const distort = (x, y) => {
      const dx = (x - c[0]) / R, dy = (y - c[1]) / R, ru = Math.hypot(dx, dy);
      if (Math.abs(lam) < 1e-12 || ru < 1e-12) return [x, y];
      const disc = 1 - 4 * lam * ru * ru;
      const rd = disc <= 0 ? 1 / (2 * lam * ru) : (1 - Math.sqrt(disc)) / (2 * lam * ru);
      const k = rd / ru;
      return [c[0] + dx * R * k, c[1] + dy * R * k];
    };
    const corrAt = corr ? (u, v) => meshAt(corr, u, v) : () => [0, 0];
    const fwd0 = (u, v) => { const p = applyH(H, u, v); return distort(p[0], p[1]); };
    const inv0 = (x, y) => { const p = undistort(x, y); return applyH(Hi, p[0], p[1]); };
    return {
      H, Hi, lam, c, R, corr,
      forward: (u, v) => { const p = fwd0(u, v), d = corrAt(u, v); return [p[0] + d[0], p[1] + d[1]]; },
      inverse: (x, y) => {
        let q = inv0(x, y);
        if (corr) for (let k = 0; k < 3; k++) { const d = corrAt(q[0], q[1]); q = inv0(x - d[0], y - d[1]); }
        return q;
      },
      undistort, fwd0,
    };
  }

  function fitWithLambda(matches, lam, c, R) {
    const m0 = makeModel([[1, 0, 0], [0, 1, 0], [0, 0, 1]], lam, c, R);
    const H = homography(matches.map(m => [m.u, m.v]), matches.map(m => m0.undistort(m.x, m.y)));
    if (!H) return null;
    const model = makeModel(H, lam, c, R);
    let se = 0;
    for (const m of matches) { const p = model.forward(m.u, m.v); se += (p[0] - m.x) ** 2 + (p[1] - m.y) ** 2; }
    model.rmsPx = Math.sqrt(se / matches.length);
    return model;
  }

  function fitModel(matches, c, R, withLens) {
    // RANSAC on the plain homography to drop mislabelled points
    let seed = 987654321;
    const rnd = () => ((seed = (seed * 1664525 + 1013904223) >>> 0) / 4294967296);
    const spacing = medianSpacing(matches);
    const thr = Math.max(1.5, 0.2 * spacing);
    let bestIn = null;
    for (let it = 0; it < 250; it++) {
      const s = [0, 0, 0, 0].map(() => matches[Math.floor(rnd() * matches.length)]);
      const H = homography(s.map(m => [m.u, m.v]), s.map(m => [m.x, m.y]));
      if (!H) continue;
      const inl = matches.filter(m => { const p = applyH(H, m.u, m.v); return Math.hypot(p[0] - m.x, p[1] - m.y) < thr * 3; });
      if (!bestIn || inl.length > bestIn.length) bestIn = inl;
    }
    if (!bestIn || bestIn.length < 12) return null;
    let inl = bestIn, model = null;
    for (let round = 0; round < 3; round++) {
      let lam = 0;
      if (withLens) {                    // golden-section search for λ
        let a = -0.6, b = 0.4;
        const g = (Math.sqrt(5) - 1) / 2, f = l => { const m = fitWithLambda(inl, l, c, R); return m ? m.rmsPx : 1e9; };
        let x1 = b - g * (b - a), x2 = a + g * (b - a), f1 = f(x1), f2 = f(x2);
        for (let k = 0; k < 30; k++) {
          if (f1 < f2) { b = x2; x2 = x1; f2 = f1; x1 = b - g * (b - a); f1 = f(x1); }
          else { a = x1; x1 = x2; f1 = f2; x2 = a + g * (b - a); f2 = f(x2); }
        }
        lam = (a + b) / 2;
        if (f(0) <= f(lam) * 1.02) lam = 0;   // only keep the lens term when it actually helps
      }
      model = fitWithLambda(inl, lam, c, R);
      if (!model) return null;
      const t = Math.max(1.0, 3 * model.rmsPx, 0.08 * spacing);
      const next = matches.filter(m => { const p = model.forward(m.u, m.v); return Math.hypot(p[0] - m.x, p[1] - m.y) < Math.min(t, 0.25 * spacing); });
      if (next.length < 12) break;
      inl = next;
    }
    model = fitWithLambda(inl, model.lam, c, R);
    model.inliers = inl;
    return model;
  }

  // ---------- smooth residual correction (bilinear mesh in cm space, image-px offsets) ----------
  function meshAt(m, u, v) {
    let fu = (u - m.u0) / m.step, fv = (v - m.v0) / m.step;
    fu = Math.max(0, Math.min(m.nu - 1.000001, fu)); fv = Math.max(0, Math.min(m.nv - 1.000001, fv));
    const i = Math.floor(fu), j = Math.floor(fv), a = fu - i, b = fv - j, k = j * m.nu + i;
    const w00 = (1 - a) * (1 - b), w10 = a * (1 - b), w01 = (1 - a) * b, w11 = a * b;
    return [m.dx[k] * w00 + m.dx[k + 1] * w10 + m.dx[k + m.nu] * w01 + m.dx[k + m.nu + 1] * w11,
            m.dy[k] * w00 + m.dy[k + 1] * w10 + m.dy[k + m.nu] * w01 + m.dy[k + m.nu + 1] * w11];
  }
  // Banded SPD solve (Cholesky). L holds the lower band: (i, j) with 0 ≤ i−j ≤ bw at i*(bw+1) + (i−j).
  function bandSolve(L, n, bw, rhs) {
    const W1 = bw + 1, at = (i, j) => i * W1 + (i - j);
    for (let j = 0; j < n; j++) {
      let d = L[at(j, j)];
      for (let k = Math.max(0, j - bw); k < j; k++) d -= L[at(j, k)] ** 2;
      if (d <= 0) return null;
      d = Math.sqrt(d); L[at(j, j)] = d;
      for (let i = j + 1; i <= Math.min(n - 1, j + bw); i++) {
        let s = L[at(i, j)];
        for (let k = Math.max(0, i - bw); k < j; k++) s -= L[at(i, k)] * L[at(j, k)];
        L[at(i, j)] = s / d;
      }
    }
    return rhs.map(b => {
      const y = Float64Array.from(b);
      for (let i = 0; i < n; i++) { for (let k = Math.max(0, i - bw); k < i; k++) y[i] -= L[at(i, k)] * y[k]; y[i] /= L[at(i, i)]; }
      for (let i = n - 1; i >= 0; i--) { for (let k = i + 1; k <= Math.min(n - 1, i + bw); k++) y[i] -= L[at(k, i)] * y[k]; y[i] /= L[at(i, i)]; }
      return y;
    });
  }
  function fitMesh(base, pts, step, alpha = 3) {
    let u0 = Infinity, u1 = -Infinity, v0 = Infinity, v1 = -Infinity;
    for (const m of pts) { u0 = Math.min(u0, m.u); u1 = Math.max(u1, m.u); v0 = Math.min(v0, m.v); v1 = Math.max(v1, m.v); }
    u0 -= step; v0 -= step;
    const nu = Math.ceil((u1 - u0) / step) + 2, nv = Math.ceil((v1 - v0) / step) + 2, n = nu * nv, bw = nu + 1;
    const L = new Float64Array(n * (bw + 1)), bx = new Float64Array(n), by = new Float64Array(n);
    const add = (i, j, w) => { if (i < j) [i, j] = [j, i]; L[i * (bw + 1) + (i - j)] += w; };
    for (const m of pts) {
      const p = base.fwd0(m.u, m.v), rx = m.x - p[0], ry = m.y - p[1];
      const fu = (m.u - u0) / step, fv = (m.v - v0) / step, i = Math.floor(fu), j = Math.floor(fv), a = fu - i, b = fv - j;
      const ks = [j * nu + i, j * nu + i + 1, (j + 1) * nu + i, (j + 1) * nu + i + 1];
      const ws = [(1 - a) * (1 - b), a * (1 - b), (1 - a) * b, a * b];
      for (let x = 0; x < 4; x++) {
        bx[ks[x]] += ws[x] * rx; by[ks[x]] += ws[x] * ry;
        for (let y = 0; y <= x; y++) add(ks[x], ks[y], x === y ? ws[x] * ws[x] : ws[x] * ws[y]);
      }
    }
    const beta = 1e-3;                     // alpha: membrane smoothness between neighbouring nodes; beta: weak pull to zero
    for (let j = 0; j < nv; j++) for (let i = 0; i < nu; i++) {
      const k = j * nu + i; add(k, k, beta);
      for (const kk of [i + 1 < nu ? k + 1 : -1, j + 1 < nv ? k + nu : -1]) {
        if (kk < 0) continue;
        add(k, k, alpha); add(kk, kk, alpha); add(kk, k, -alpha);
      }
    }
    const sol = bandSolve(L, n, bw, [bx, by]);
    if (!sol) return null;
    return { u0, v0, step, nu, nv, dx: sol[0], dy: sol[1] };
  }
  const rmsOf = (model, pts) => {
    let se = 0; for (const m of pts) { const p = model.forward(m.u, m.v); se += (p[0] - m.x) ** 2 + (p[1] - m.y) ** 2; }
    return Math.sqrt(se / pts.length);
  };
  // Pick the mesh spacing that best predicts held-out crossings (2-fold); none if no spacing helps.
  function withMesh(model, inl) {
    if (inl.length < 80) return model;
    const A = inl.filter((m, k) => k % 2 === 0), B = inl.filter((m, k) => k % 2 === 1);
    const score = (step, alpha) => {
      let e = 0;
      for (const [fitSet, testSet] of [[A, B], [B, A]]) {
        const mesh = step ? fitMesh(model, fitSet, step, alpha) : null;
        if (step && !mesh) return Infinity;
        e += rmsOf(step ? makeModel(model.H, model.lam, model.c, model.R, mesh) : model, testSet);
      }
      return e;
    };
    const base = score(0);
    let best = null, bestE = base;
    for (const step of [6, 4, 3, 2]) for (const alpha of [3, 0.5]) {
      const e = score(step, alpha); if (e < bestE) { bestE = e; best = [step, alpha]; }
    }
    if (!best || bestE > 0.95 * base) return model;
    const out = makeModel(model.H, model.lam, model.c, model.R, fitMesh(model, inl, ...best));
    out.rmsPx = rmsOf(out, inl); out.inliers = inl; out.meshGain = base / bestE; out.meshStep = best[0];
    return out;
  }

  // The sweep validates every crossing on the image itself, so keep crossings the global model
  // can't explain if a smooth local mesh can; only trim points that disagree with their neighbourhood.
  // A real crossing agrees with its neighbours: predicted from each found neighbour (neighbour +
  // model step) it lands within ~1 mm. False crossings (radial "starburst" lines, text, paper
  // marks) don't; nor do isolated singles, which can't be verified.
  function neighbourConsistent(model, pts, spacing) {
    const byK = new Map(pts.map(m => [m.u + ',' + m.v, m])), tol = 0.1 * spacing;
    return pts.filter(m => {
      const pm = model.forward(m.u, m.v);
      let agree = 0;
      for (const [du, dv] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
        const n = byK.get((m.u + du) + ',' + (m.v + dv)); if (!n) continue;
        const pn = model.forward(n.u, n.v);
        if (Math.hypot(m.x - (n.x + pm[0] - pn[0]), m.y - (n.y + pm[1] - pn[1])) < tol) agree++;
      }
      return agree >= 2;
    });
  }

  function finalFit(model, all) {
    const spacing = medianSpacing(all);
    const checked = neighbourConsistent(model, all, spacing);
    if (checked.length >= 30) all = checked;
    let pts = all;
    for (let k = 0; k < 2; k++) {
      const mesh = fitMesh(model, pts, 2, 0.5);
      if (!mesh) break;
      const mm = makeModel(model.H, model.lam, model.c, model.R, mesh);
      const res = all.map(m => { const p = mm.forward(m.u, m.v); return Math.hypot(p[0] - m.x, p[1] - m.y); });
      const med = [...res].sort((a, b) => a - b)[res.length >> 1];
      // a mislabelled crossing is off by a whole cell; real ones stay within a fraction of it
      const keep = all.filter((m, i) => res[i] < Math.max(0.35 * spacing, 5 * 1.4826 * med));
      if (keep.length < 12) break;
      pts = keep;
    }
    const out = withMesh(model, pts);
    if (!out.corr) return model;             // mesh didn't help on held-out points: global model only
    return out;
  }

  function medianSpacing(matches) {
    const byKey = new Map(matches.map(m => [m.u + ',' + m.v, m]));
    const d = [];
    for (const m of matches) {
      const n = byKey.get((m.u + 1) + ',' + m.v);
      if (n) d.push(Math.hypot(n.x - m.x, n.y - m.y));
      const n2 = byKey.get(m.u + ',' + (m.v + 1));
      if (n2) d.push(Math.hypot(n2.x - m.x, n2.y - m.y));
    }
    d.sort((a, b) => a - b);
    return d.length ? d[d.length >> 1] : 10;
  }

  // ---------- image helpers ----------
  // intensity-weighted centroid of the line image around (x,y); converges onto line crossings
  function centroid(t, W, H, x, y, r, iters = 6) {
    for (let it = 0; it < iters; it++) {
      const xi = Math.round(x), yi = Math.round(y);
      let sx = 0, sy = 0, sw = 0;
      for (let dy = -r; dy <= r; dy++) {
        const py = yi + dy; if (py < 0 || py >= H) continue;
        for (let dx = -r; dx <= r; dx++) {
          const px = xi + dx; if (px < 0 || px >= W) continue;
          let w = t[py * W + px]; if (!w) continue;
          w *= w; sx += w * px; sy += w * py; sw += w;
        }
      }
      if (!sw) return null;
      const nx = sx / sw, ny = sy / sw;
      const moved = Math.hypot(nx - x, ny - y);
      x = nx; y = ny;
      if (moved < 0.02) break;
    }
    return [x, y];
  }
  const sampleMax = (t, W, H, x, y) => {     // max of 3×3 around (x,y)
    let m = 0; const xi = Math.round(x), yi = Math.round(y);
    for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
      const px = xi + dx, py = yi + dy;
      if (px >= 0 && py >= 0 && px < W && py < H) m = Math.max(m, t[py * W + px]);
    }
    return m;
  };

  const bil = (t, W, H, x, y) => {
    const x0 = Math.floor(x), y0 = Math.floor(y);
    if (x0 < 0 || y0 < 0 || x0 >= W - 1 || y0 >= H - 1) return 0;
    const fx = x - x0, fy = y - y0, i = y0 * W + x0;
    return t[i] * (1 - fx) * (1 - fy) + t[i + 1] * fx * (1 - fy) + t[i + W] * (1 - fx) * fy + t[i + W + 1] * fx * fy;
  };
  // Precise crossing: locate each of the two lines on its arms (away from the crossing blob)
  // with sub-pixel profile centroids, fit both lines, intersect them.
  function refineCrossing(t, W, H, q, e1, e2, maxSlope = 0.15) {
    const L1 = Math.hypot(...e1), L2 = Math.hypot(...e2);
    const d1 = [e1[0] / L1, e1[1] / L1], d2 = [e2[0] / L2, e2[1] / L2];
    const n1 = [-d1[1], d1[0]], n2 = [-d2[1], d2[0]];
    const fit = (d, n, along, across) => {
      const S = [], O = [], h = 0.3 * across, step = 0.25;
      for (let f = 0.28; f <= 0.461; f += 0.03) for (const sg of [-1, 1]) {
        const sa = sg * f * along, bx = q[0] + d[0] * sa, by = q[1] + d[1] * sa;
        const vals = [];
        for (let o = -h; o <= h + 1e-9; o += step) vals.push([o, bil(t, W, H, bx + n[0] * o, by + n[1] * o)]);
        // line centre = centroid of the line's core: subtract a linear background (glare, shading)
        // through the profile ends, then use only samples above half of the peak
        const m = vals.length, e = Math.max(2, Math.round(m * 0.12));
        let l0 = 0, r0 = 0, lo = 0, ro = 0;
        for (let k = 0; k < e; k++) { l0 += vals[k][1]; lo += vals[k][0]; r0 += vals[m - 1 - k][1]; ro += vals[m - 1 - k][0]; }
        l0 /= e; lo /= e; r0 /= e; ro /= e;
        const w = vals.map(([o, v]) => v - (l0 + (r0 - l0) * (o - lo) / (ro - lo)));
        let pk = 0; for (let k = 1; k < m; k++) if (w[k] > w[pk]) pk = k;
        if (!(w[pk] > 0)) continue;
        const half = 0.5 * w[pk];
        let a0 = pk, a1 = pk;
        while (a0 > 0 && w[a0 - 1] > half) a0--;
        while (a1 < m - 1 && w[a1 + 1] > half) a1++;
        if (a0 === 0 || a1 === m - 1) continue;          // line not contained in the window
        let sw = 0, so = 0;
        for (let k = a0; k <= a1; k++) { const ww = w[k] - half; sw += ww; so += ww * vals[k][0]; }
        if (sw > 0) { S.push(sa); O.push(so / sw); }
      }
      if (S.length < 6) return null;
      const ls = idx => {
        let ms = 0, mo = 0; for (const k of idx) { ms += S[k]; mo += O[k]; }
        ms /= idx.length; mo /= idx.length;
        let sso = 0, sss = 0; for (const k of idx) { sso += (S[k] - ms) * (O[k] - mo); sss += (S[k] - ms) ** 2; }
        const bb = sss ? sso / sss : 0;
        return { a: mo - bb * ms, b: bb };
      };
      let L = ls(S.map((_, k) => k));
      // one robust pass: drop samples disturbed by printed text or ticks next to the line
      const res = S.map((sv, k) => Math.abs(O[k] - (L.a + L.b * sv)));
      const med = [...res].sort((x, y) => x - y)[res.length >> 1];
      const keep = S.map((_, k) => k).filter(k => res[k] <= Math.max(0.15, 3 * med));
      if (keep.length >= 5 && keep.length < S.length) L = ls(keep);
      return L;
    };
    const A = fit(d1, n1, L1, L2), B = fit(d2, n2, L2, L1);
    if (!A || !B) return null;
    // both lines must run along the grid directions: rejects crossings with the mat's diagonal
    // angle lines and the radial "starburst" lines (slope b = tan of the angle off the grid axis)
    if (Math.abs(A.b) > maxSlope || Math.abs(B.b) > maxSlope) { rej('angle'); return null; }
    // q + s·d1 + (a1 + b1 s)·n1 = q + t·d2 + (a2 + b2 t)·n2
    const c1 = [d1[0] + A.b * n1[0], d1[1] + A.b * n1[1]], c2 = [-(d2[0] + B.b * n2[0]), -(d2[1] + B.b * n2[1])];
    const rhs = [B.a * n2[0] - A.a * n1[0], B.a * n2[1] - A.a * n1[1]];
    const det = c1[0] * c2[1] - c1[1] * c2[0];
    if (Math.abs(det) < 1e-9) return null;
    const ss = (rhs[0] * c2[1] - rhs[1] * c2[0]) / det;
    const x = q[0] + ss * d1[0] + (A.a + A.b * ss) * n1[0], y = q[1] + ss * d1[1] + (A.a + A.b * ss) * n1[1];
    if (Math.hypot(x - q[0], y - q[1]) > 0.2 * Math.min(L1, L2)) return null;
    return [x, y];
  }

  // ---------- lattice growing from detected corner candidates ----------
  // Lock onto the grid without any model: find a seed whose 3×3 neighbourhood is 9 validated
  // crossings, then grow by single validated steps (each node predicts its neighbour from its own
  // measured steps). Corner candidates are only used to propose seeds.
  function buildLattice(pts, t, W, H, centre, thrMin) {
    const cell = 24, cols = Math.ceil(W / cell), rows = Math.ceil(H / cell);
    const buckets = new Map();
    pts.forEach((p, i) => {
      const k = Math.floor(p[1] / cell) * cols + Math.floor(p[0] / cell);
      if (!buckets.has(k)) buckets.set(k, []); buckets.get(k).push(i);
    });
    function neighbours(p, rmin, rmax) {
      const out = [];
      const c0 = Math.floor((p[0] - rmax) / cell), c1 = Math.floor((p[0] + rmax) / cell), r0 = Math.floor((p[1] - rmax) / cell), r1 = Math.floor((p[1] + rmax) / cell);
      for (let rr = Math.max(0, r0); rr <= Math.min(rows - 1, r1); rr++)
        for (let cc = Math.max(0, c0); cc <= Math.min(cols - 1, c1); cc++) {
          const b = buckets.get(rr * cols + cc); if (!b) continue;
          for (const i of b) { const d = dist(p, pts[i]); if (d >= rmin && d <= rmax) out.push([i, d]); }
        }
      return out.sort((a, b) => a[1] - b[1]);
    }
    // seed: candidates nearest the mat centre first; accept the first clean 3×3 block
    const order = pts.map((p, i) => i).sort((a, b) => dist(pts[a], centre) - dist(pts[b], centre));
    function nextSeed(start, skip) {
    for (let oi = start; oi < Math.min(order.length, 3000); oi++) {
      const i = order[oi];
      if (skip(pts[i])) continue;
      const p0 = pts[i], nb = neighbours(p0, 5, 90);
      if (nb.length < 2) continue;
      const q1 = pts[nb[0][0]], a = [q1[0] - p0[0], q1[1] - p0[1]], la = Math.hypot(...a);
      let b = null;
      for (const [j, d] of nb) {
        if (d > 1.5 * la || d < 0.66 * la) continue;
        const v = [pts[j][0] - p0[0], pts[j][1] - p0[1]];
        if (Math.abs((v[0] * a[0] + v[1] * a[1]) / (d * la)) < 0.35) { b = v; break; }
      }
      if (!b) continue;
      // candidate-based directions are rough here: allow ~20° until the block is measured
      const c = findCrossing(t, W, H, p0, a, b, thrMin, 0.35);
      if (!c) continue;
      const block = new Map();
      let ok = true;
      for (let j = -1; j <= 1 && ok; j++) for (let k = -1; k <= 1 && ok; k++) {
        const q = (j || k) ? findCrossing(t, W, H, [c[0] + j * a[0] + k * b[0], c[1] + j * a[1] + k * b[1]], a, b, thrMin, 0.35) : c;
        if (!q) ok = false; else block.set(j + ',' + k, q);
      }
      if (!ok) continue;
      // the block must be a regular lattice (rules out radial-line junk)
      const s1 = [block.get('1,0')[0] - c[0], block.get('1,0')[1] - c[1]], s2 = [block.get('0,1')[0] - c[0], block.get('0,1')[1] - c[1]];
      const d11 = block.get('1,1'), pred = [c[0] + s1[0] + s2[0], c[1] + s1[1] + s2[1]];
      if (Math.hypot(d11[0] - pred[0], d11[1] - pred[1]) > 0.15 * Math.min(Math.hypot(...s1), Math.hypot(...s2))) continue;
      return { c, a: s1, b: s2, oi };
    }
    return null;
    }
    // grow from a seed; try a few seeds (clutter like the radial "starburst" can box one in) and
    // keep the largest lattice
    function growFrom(seed) {
      let e1 = seed.a, e2 = seed.b;
      // orient: e1 ≈ image +x, e2 ≈ image +y
      if (Math.abs(e1[0]) < Math.abs(e2[0])) [e1, e2] = [e2, e1];
      if (e1[0] < 0) e1 = [-e1[0], -e1[1]];
      if (e2[1] < 0) e2 = [-e2[0], -e2[1]];

      const nodes = new Map(), key = (i, j) => i + ',' + j, queue = [];
      const hash = new Map(), hk = (x, y) => Math.floor(y / 8) * 100000 + Math.floor(x / 8);
      const taken = (x, y, r) => {
        for (let gy = Math.floor((y - r) / 8); gy <= Math.floor((y + r) / 8); gy++)
          for (let gx = Math.floor((x - r) / 8); gx <= Math.floor((x + r) / 8); gx++) {
            const b = hash.get(gy * 100000 + gx); if (!b) continue;
            for (const n of b) if (Math.hypot(n.x - x, n.y - y) < r) return true;
          }
        return false;
      };
      const add = (i, j, q, a, b) => {
        const n = { i, j, x: q[0], y: q[1], e1: a, e2: b }; nodes.set(key(i, j), n); queue.push(n);
        const k = hk(q[0], q[1]); if (!hash.has(k)) hash.set(k, []); hash.get(k).push(n);
      };
      add(0, 0, seed.c, e1, e2);
      while (queue.length && nodes.size < 30000) {
        const n = queue.shift();
        for (const [di, dj] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
          const ti = n.i + di, tj = n.j + dj;
          if (nodes.has(key(ti, tj))) continue;
          const back = nodes.get(key(n.i - di, n.j - dj));
          const base = di ? [n.e1[0] * di, n.e1[1] * di] : [n.e2[0] * dj, n.e2[1] * dj];
          let e = back ? [n.x - back.x, n.y - back.y] : base;
          if (Math.abs(Math.hypot(...e) / Math.hypot(...base) - 1) > 0.3) e = base;
          const ae1 = di ? [e[0] * di, e[1] * di] : n.e1, ae2 = dj ? [e[0] * dj, e[1] * dj] : n.e2;
          const q = findCrossing(t, W, H, [n.x + e[0], n.y + e[1]], ae1, ae2, thrMin);
          if (!q || taken(q[0], q[1], 0.3 * Math.min(Math.hypot(...ae1), Math.hypot(...ae2)))) continue;
          const step = [q[0] - n.x, q[1] - n.y];
          if (Math.hypot(step[0] - e[0], step[1] - e[1]) > 0.2 * Math.hypot(...e)) continue;   // irregular: not the grid
          add(ti, tj, q, di ? [step[0] * di, step[1] * di] : n.e1, dj ? [step[0] * dj, step[1] * dj] : n.e2);
        }
      }
      return [...nodes.values()].map(n => ({ u: n.i, v: n.j, x: n.x, y: n.y }));
    }
    let best = null, oi = 0;
    const tried = [];
    for (let attempt = 0; attempt < 8; attempt++) {
      const seed = nextSeed(oi, p => tried.some(l => l.some(m => Math.hypot(m.x - p[0], m.y - p[1]) < 30)));
      if (!seed) break;
      oi = seed.oi + 1;
      const lat = growFrom(seed);
      tried.push(lat);
      if (!best || lat.length > best.length) best = lat;
      if (best.length >= 300) break;
    }
    return best;
  }

  // ---------- crossing search guided by the fitted model ----------
  const stats = {}, debug = {};
  const rej = k => { stats[k] = (stats[k] || 0) + 1; };
  // Look for a real grid crossing near p (e1/e2 = one-cell steps). Thresholds are local, so glare or
  // shadow across the mat doesn't matter; thrMin only rejects near-empty areas.
  function findCrossing(t, W, H, p, e1, e2, thrMin, maxSlope) {
    if (!(p[0] > 2 && p[1] > 2 && p[0] < W - 3 && p[1] < H - 3)) return null;
    const s = Math.min(Math.hypot(...e1), Math.hypot(...e2));
    if (!(s > 4)) return null;
    const q = centroid(t, W, H, p[0], p[1], Math.max(2, Math.round(0.3 * s)));
    if (!q || Math.hypot(q[0] - p[0], q[1] - p[1]) > 0.25 * s) { rej('centroid'); return null; }
    const peak = sampleMax(t, W, H, q[0], q[1]);
    if (peak < thrMin) { rej('faint'); return null; }
    const thr = 0.4 * peak;
    // bright along all 4 arms, dark in the middle of the 4 neighbouring squares
    const arm = (dx, dy) => sampleMax(t, W, H, q[0] + dx, q[1] + dy) >= thr;
    if (!(arm(0.4 * e1[0], 0.4 * e1[1]) && arm(-0.4 * e1[0], -0.4 * e1[1]) && arm(0.4 * e2[0], 0.4 * e2[1]) && arm(-0.4 * e2[0], -0.4 * e2[1]))) { rej('arms'); return null; }
    let dark = 0;
    for (const [a, b] of [[0.5, 0.5], [-0.5, 0.5], [0.5, -0.5], [-0.5, -0.5]]) {
      const x = Math.round(q[0] + a * e1[0] + b * e2[0]), y = Math.round(q[1] + a * e1[1] + b * e2[1]);
      if (x >= 0 && y >= 0 && x < W && y < H && t[y * W + x] < thr) dark++;
    }
    if (dark < 3) { rej('dark'); return null; }
    const fine = refineCrossing(t, W, H, q, e1, e2, maxSlope);
    if (!fine) { rej('refine'); return null; }
    if (Math.hypot(fine[0] - p[0], fine[1] - p[1]) > 0.25 * s) { rej('drift'); return null; }
    return fine;
  }

  // Label-consistent crossing search.
  // Crossings are only labelled by 1-cell steps from an already-found neighbour (validated within
  // 0.25 cell), so every connected component is internally consistent. Component 0 grows from the
  // trusted seeds. Areas it can't reach (cut off by paper, weights, tape) are seeded from the model's
  // predictions and grown as separate components; when one runs into an already-found crossing we
  // learn its exact label offset, otherwise the offset comes from the model (dropped if ambiguous).
  function sweep(model, t, W, H, uvBox, thrMin, seeds) {
    for (const k in stats) delete stats[k];
    const cell = 8, hash = new Map();
    const hk = (x, y) => Math.floor(y / cell) * 100000 + Math.floor(x / cell);
    const nearPos = (x, y, r) => {
      let best = null, bd = r;
      for (let gy = Math.floor((y - r) / cell); gy <= Math.floor((y + r) / cell); gy++)
        for (let gx = Math.floor((x - r) / cell); gx <= Math.floor((x + r) / cell); gx++) {
          const b = hash.get(gy * 100000 + gx); if (!b) continue;
          for (const n of b) { const d = Math.hypot(n.x - x, n.y - y); if (d < bd) { bd = d; best = n; } }
        }
      return best;
    };
    const steps = (u, v, p) => {
      const pu = model.forward(u + 1, v), pv = model.forward(u, v + 1);
      return [[pu[0] - p[0], pu[1] - p[1]], [pv[0] - p[0], pv[1] - p[1]]];
    };
    const comps = [];
    const newComp = () => { const c = { id: comps.length, nodes: new Map(), votes: new Map() }; comps.push(c); return c; };
    const addNode = (c, u, v, x, y) => {
      const n = { u, v, x, y, c: c.id }; c.nodes.set(u + ',' + v, n);
      const k = hk(x, y); if (!hash.has(k)) hash.set(k, []); hash.get(k).push(n);
      return n;
    };
    // growth is bounded by the green mat itself (no crossings elsewhere); this is only a safety cap
    const inBox = (u, v) => Math.abs(u) < 400 && Math.abs(v) < 400;
    function grow(c, start) {
      const queue = [...start];
      while (queue.length) {
        const m = queue.shift();
        const pm = model.forward(m.u, m.v);
        for (const [du, dv] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
          const u = m.u + du, v = m.v + dv;
          if (c.nodes.has(u + ',' + v) || !inBox(u, v)) continue;
          const pt = model.forward(u, v);
          const back = c.nodes.get((m.u - du) + ',' + (m.v - dv));
          const p = back ? [2 * m.x - back.x, 2 * m.y - back.y] : [m.x + pt[0] - pm[0], m.y + pt[1] - pm[1]];
          let [e1, e2] = steps(u, v, pt);
          if (back) { const st = [(m.x - back.x) * (du || dv), (m.y - back.y) * (du || dv)]; if (du) e1 = st; else e2 = st; }
          const q = findCrossing(t, W, H, p, e1, e2, thrMin);
          if (!q) continue;
          const s = Math.min(Math.hypot(...e1), Math.hypot(...e2));
          if (back && Math.hypot(q[0] - p[0], q[1] - p[1]) > 0.2 * s) continue;   // irregular step: not the grid
          const other = nearPos(q[0], q[1], 0.3 * s);
          if (other) {                                  // physical crossing already labelled
            if (other.c !== c.id) {
              const vk = other.c + ':' + (other.u - u) + ',' + (other.v - v);
              c.votes.set(vk, (c.votes.get(vk) || 0) + 1);
            }
            continue;
          }
          queue.push(addNode(c, u, v, q[0], q[1]));
        }
      }
    }
    // component 0 from the trusted seeds
    const c0 = newComp(), start = [];
    for (const m of seeds) if (!c0.nodes.has(m.u + ',' + m.v) && !nearPos(m.x, m.y, 2)) start.push(addNode(c0, m.u, m.v, m.x, m.y));
    grow(c0, start);
    // other areas, seeded where the model's prediction lands right on a crossing
    for (let v = uvBox.v0; v <= uvBox.v1; v++) for (let u = uvBox.u0; u <= uvBox.u1; u++) {
      const p = model.forward(u, v), [e1, e2] = steps(u, v, p);
      const s = Math.min(Math.hypot(...e1), Math.hypot(...e2));
      if (!(s > 4) || nearPos(p[0], p[1], 0.6 * s)) continue;
      const q = findCrossing(t, W, H, p, e1, e2, thrMin);
      if (!q || nearPos(q[0], q[1], 0.3 * s)) continue;
      const c = newComp();
      grow(c, [addNode(c, u, v, q[0], q[1])]);
    }
    // 1) merge components that ran into each other (exact offsets), strongest links first
    const parent = comps.map((c, i) => i), off = comps.map(() => [0, 0]);   // label in parent = label + off
    const find = i => {
      if (parent[i] === i) return i;
      const p = parent[i], r = find(p);
      off[i] = [off[i][0] + off[p][0], off[i][1] + off[p][1]]; parent[i] = r;
      return r;
    };
    const links = [];
    for (const c of comps) {
      const best = new Map();
      for (const [vk, n] of c.votes) { const [oc] = vk.split(':'); if (!best.has(oc) || best.get(oc)[1] < n) best.set(oc, [vk, n]); }
      for (const [oc, [vk, n]] of best) { const d = vk.split(':')[1].split(',').map(Number); links.push({ a: c.id, b: +oc, d, n }); }
    }
    links.sort((x, y) => y.n - x.n);
    for (const { a, b, d } of links) {
      const ra = find(a), rb = find(b);
      if (ra === rb) continue;
      off[ra] = [-off[a][0] + d[0] + off[b][0], -off[a][1] + d[1] + off[b][1]]; parent[ra] = rb;
    }
    const groups = new Map();
    for (const c of comps) {
      const r = find(c.id), o = off[c.id];
      if (!groups.has(r)) groups.set(r, []);
      for (const n of c.nodes.values()) groups.get(r).push({ u: n.u + o[0], v: n.v + o[1], x: n.x, y: n.y });
    }
    // 2) the largest merged group is the reference frame; 3) place the others with a model fitted
    //    to everything placed so far, refitting as groups are added; ambiguous groups are dropped
    const gl = [...groups.values()].filter(g => g.length >= 6).sort((x, y) => y.length - x.length);
    if (!gl.length) return [];
    let placed = gl[0], pending = gl.slice(1);
    for (let progress = true; progress && pending.length;) {
      progress = false;
      const m = fitModel(placed, model.c, model.R, true) || model;
      const still = [];
      for (const g of pending) {
        const du = [], dv = [];
        for (const n of g) { const q = m.inverse(n.x, n.y); du.push(q[0] - n.u); dv.push(q[1] - n.v); }
        du.sort((a, b) => a - b); dv.sort((a, b) => a - b);
        const mu = du[du.length >> 1], mv = dv[dv.length >> 1];
        if (Math.abs(mu - Math.round(mu)) > 0.3 || Math.abs(mv - Math.round(mv)) > 0.3) { still.push(g); continue; }
        const o = [Math.round(mu), Math.round(mv)];
        placed = placed.concat(g.map(n => ({ u: n.u + o[0], v: n.v + o[1], x: n.x, y: n.y })));
        progress = true;
      }
      pending = still;
    }
    if (pending.length) rej('ambiguous-region');
    const out = placed;
    // a label claimed twice (shouldn't happen after resolution): keep the first (component 0 first)
    const seen = new Set();
    const res = out.filter(m => { const k = m.u + ',' + m.v; if (seen.has(k)) return false; seen.add(k); return true; });
    debug.components = comps.length;
    debug.groups = gl.map(g => g.length); debug.dropped = pending.map(g => g.length);
    return res;
  }

  // ---------- main entry ----------
  // bitmap → { model (analysis coords), sc/scY (analysis px per photo px, per axis), W, H, matches, rmsMm, spacingPx, visible:{u0,u1,v0,v1} }
  // Very close shots draw the grid with huge squares (tens of px per cm), which the line filters
  // aren't sized for: if locking on fails, retry on a smaller analysis image.
  function detect(bitmap, maskFn, maxDim = 1600) {
    const r = detectAt(bitmap, maskFn, maxDim);
    if (r.ok || maxDim <= 700) return r;
    const r2 = detectAt(bitmap, maskFn, Math.round(maxDim / 2));
    return r2.ok ? r2 : r;
  }
  function detectAt(bitmap, maskFn, maxDim) {
    const scale = Math.min(1, maxDim / Math.max(bitmap.width, bitmap.height));
    const c = document.createElement('canvas');
    c.width = Math.round(bitmap.width * scale); c.height = Math.round(bitmap.height * scale);
    c.getContext('2d').drawImage(bitmap, 0, 0, c.width, c.height);
    // per-axis scales: width and height are rounded separately, so they can differ slightly
    const sc = c.width / bitmap.width, scY = c.height / bitmap.height, W = c.width, H = c.height;
    const src = cv.imread(c);
    const mask = maskFn(src);
    const gray = new cv.Mat(), th = new cv.Mat(), tb = new cv.Mat(), masked = cv.Mat.zeros(H, W, cv.CV_8UC1);
    cv.cvtColor(src, gray, cv.COLOR_RGBA2GRAY);
    const k = cv.getStructuringElement(cv.MORPH_RECT, new cv.Size(11, 11));
    cv.morphologyEx(gray, th, cv.MORPH_TOPHAT, k);
    th.copyTo(masked, mask);
    cv.GaussianBlur(masked, tb, new cv.Size(3, 3), 0);
    const t = new Float32Array(tb.data);
    // corner candidates
    const corners = new cv.Mat();
    cv.goodFeaturesToTrack(tb, corners, 40000, 0.01, 3, mask, 3, false, 0.04);
    let pts = [];
    for (let i = 0; i < corners.rows; i++) {
      const q = centroid(t, W, H, corners.data32F[2 * i], corners.data32F[2 * i + 1], 3, 4);
      if (q) pts.push(q);
    }
    // mat region outline (for centre + visible extent)
    const contours = new cv.MatVector(), hier = new cv.Mat();
    cv.findContours(mask, contours, hier, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_NONE);
    const regions = [];
    for (let i = 0; i < contours.size(); i++) {
      const cn = contours.get(i);
      if (cv.contourArea(cn) < 0.002 * W * H) continue;
      const poly = [];
      for (let j = 0; j < cn.data32S.length; j += 2) poly.push([cn.data32S[j], cn.data32S[j + 1]]);
      regions.push(poly);
    }
    let outline = regions.flat();
    [src, mask, gray, th, tb, masked, k, corners, contours, hier].forEach(m => m.delete());

    const fail = reason => ({ ok: false, reason, sc, scY, W, H, candidates: pts.length });
    if (pts.length < 30 || !outline.length) return fail('no grid lines found on the green area');
    let mx = 0, my = 0; for (const p of outline) { mx += p[0]; my += p[1]; }
    const lattice = buildLattice(pts, t, W, H, [mx / outline.length, my / outline.length], 10);
    if (!lattice || lattice.length < 30) return fail('could not lock onto the 1 cm grid');
    // only green regions that actually carry the grid count as mat (ignores other green objects)

    const cImg = [W / 2, H / 2], R = Math.hypot(W, H) / 2;
    let model = fitModel(lattice, cImg, R, false);
    if (!model) return fail('grid fit failed');
    // thresholds for crossing validation: typical line brightness
    const vals = model.inliers.map(m => t[Math.round(m.y) * W + Math.round(m.x)]).sort((a, b) => a - b);
    const thrT = 0.25 * vals[vals.length >> 1];
    let matches = lattice; debug.rounds = [];
    for (let round = 0; round < 3; round++) {
      const box = uvBounds(model, outline, matches, 2, 40);
      const sw = sweep(model, t, W, H, box, thrT, model.inliers);
      debug.rounds.push({ n: sw.length, rmsPx: model.rmsPx, lam: model.lam });
      if (sw.length < 30) break;
      matches = sw;
      const m2 = fitModel(matches, cImg, R, true);
      if (!m2) break;
      model = m2;
    }
    model = finalFit(model, matches);
    const inl = model.inliers;
    // error in cm space
    let se = 0; const errs = [];
    for (const m of inl) { const q = model.inverse(m.x, m.y); const e = Math.hypot(q[0] - m.u, q[1] - m.v); se += e * e; errs.push(e); }
    errs.sort((a, b) => a - b);
    // only green regions that actually carry the grid count as mat (ignores other green objects)
    const onMat = regions.filter(poly => inl.filter(m => inPoly(poly, m.x, m.y)).length >= 5);
    if (onMat.length) outline = onMat.flat();
    const visible = uvBounds(model, outline, inl, 0, 15);
    // the photo's own border is not a mat edge
    const edges = matEdges(model, outline.filter(p => p[0] > 2 && p[1] > 2 && p[0] < W - 3 && p[1] < H - 3), inl);
    return {
      ok: true, model, sc, scY, W, H, matches: inl, candidates: pts.length,
      rmsMm: Math.sqrt(se / inl.length) * 10, p95Mm: errs[Math.floor(errs.length * 0.95)] * 10,
      spacingPx: medianSpacing(inl) / sc, visible, edges, lam: model.lam, outline, mesh: model.meshStep || 0,
    };
  }

  // Mat edges in cm: the mat's sides are straight, axis-aligned lines in grid coordinates, so each
  // one shows up as a sharp peak in the histogram of the green outline's u (or v) values. Weights
  // and other occluders only add scattered values. Returns the found edges (null where hidden).
  function matEdges(model, outline, matches) {
    const us = matches.map(m => m.u).sort((a, b) => a - b), vs = matches.map(m => m.v).sort((a, b) => a - b);
    const cu = us[us.length >> 1], cv = vs[vs.length >> 1];
    const pu = [], pv = [];
    for (const p of outline) { const q = model.inverse(p[0], p[1]); if (isFinite(q[0]) && isFinite(q[1])) { pu.push(q[0]); pv.push(q[1]); } }
    const spanU = us[us.length - 1] - us[0], spanV = vs[vs.length - 1] - vs[0];
    const find = (vals, other, c, side, otherSpan) => {
      const bins = new Map();
      vals.forEach(x => { if ((x - c) * side > 1) { const b = Math.round(x * 5); bins.set(b, (bins.get(b) || 0) + 1); } });
      let best = null, bs = 0;
      for (const [b] of bins) {
        const sc = (bins.get(b - 1) || 0) + bins.get(b) + (bins.get(b + 1) || 0);
        if (sc > bs || (sc === bs && (b - best) * side > 0)) { bs = sc; best = b; }
      }
      if (best === null) return null;
      const near = [], along = [];
      vals.forEach((x, i) => { if (Math.abs(x - best / 5) < 0.4) { near.push(x); along.push(other[i]); } });
      along.sort((a, b) => a - b);
      // a real edge is a long straight run: it must cover a good part of the mat's other dimension
      if (!along.length || along[along.length - 1] - along[0] < 0.25 * otherSpan) return null;
      near.sort((a, b) => a - b);
      const e = near[near.length >> 1];
      // the mat's outer edge has no grid beyond it (a straight paper edge does)
      const beyond = matches.filter(m => ((vals === pu ? m.u : m.v) - e) * side > 1).length;
      if (beyond > Math.max(5, 0.01 * matches.length)) return null;
      return e;
    };
    return {
      u0: find(pu, pv, cu, -1, spanV), u1: find(pu, pv, cu, 1, spanV),
      v0: find(pv, pu, cv, -1, spanU), v1: find(pv, pu, cv, 1, spanU),
    };
  }

  function inPoly(poly, x, y) {
    let inside = false;
    for (let i = 0, j = poly.length - 1; i < poly.length; j = i++) {
      const [xi, yi] = poly[i], [xj, yj] = poly[j];
      if ((yi > y) !== (yj > y) && x < ((xj - xi) * (y - yi)) / (yj - yi) + xi) inside = !inside;
    }
    return inside;
  }

  // integer cm box covering the visible mat (outline mapped through the model), limited to
  // a margin around the matched points so extrapolation can't blow it up
  function uvBounds(model, outline, matches, pad, lim) {
    let u0 = Infinity, u1 = -Infinity, v0 = Infinity, v1 = -Infinity;
    for (const p of outline) {
      const q = model.inverse(p[0], p[1]);
      if (!isFinite(q[0]) || !isFinite(q[1])) continue;
      u0 = Math.min(u0, q[0]); u1 = Math.max(u1, q[0]); v0 = Math.min(v0, q[1]); v1 = Math.max(v1, q[1]);
    }
    let mu0 = Infinity, mu1 = -Infinity, mv0 = Infinity, mv1 = -Infinity;
    for (const m of matches) { mu0 = Math.min(mu0, m.u); mu1 = Math.max(mu1, m.u); mv0 = Math.min(mv0, m.v); mv1 = Math.max(mv1, m.v); }
    return {
      u0: Math.floor(Math.max(u0, mu0 - lim)) - pad, u1: Math.ceil(Math.min(u1, mu1 + lim)) + pad,
      v0: Math.floor(Math.max(v0, mv0 - lim)) - pad, v1: Math.ceil(Math.min(v1, mv1 + lim)) + pad,
    };
  }

  // output pixel grid → analysis-image coordinates, written into ax/ay (Float32Array) for rows y0..y0+rows
  function fillMaps(model, r, S, y0, rows, W, ax, ay, sc, scY = sc) {
    const [[h00, h01, h02], [h10, h11, h12], [h20, h21, h22]] = model.H, lam = model.lam, [cx, cy] = model.c, R = model.R, m = model.corr;
    for (let yy = 0, i = 0; yy < rows; yy++) {
      const v = r.v0 + (y0 + yy + 0.5) / S;
      for (let x = 0; x < W; x++, i++) {
        const u = r.u0 + (x + 0.5) / S;
        const w = h20 * u + h21 * v + h22;
        let px = (h00 * u + h01 * v + h02) / w, py = (h10 * u + h11 * v + h12) / w;
        if (lam) {
          const dx = (px - cx) / R, dy = (py - cy) / R, ru = Math.sqrt(dx * dx + dy * dy);
          if (ru > 1e-12) {
            const disc = 1 - 4 * lam * ru * ru;
            const k = (disc <= 0 ? 1 / (2 * lam * ru) : (1 - Math.sqrt(disc)) / (2 * lam * ru)) / ru;
            px = cx + dx * R * k; py = cy + dy * R * k;
          }
        }
        if (m) {                                     // bilinear mesh offset, inlined for speed
          let fu = (u - m.u0) / m.step, fv = (v - m.v0) / m.step;
          fu = fu < 0 ? 0 : fu > m.nu - 1.000001 ? m.nu - 1.000001 : fu;
          fv = fv < 0 ? 0 : fv > m.nv - 1.000001 ? m.nv - 1.000001 : fv;
          const ii = fu | 0, jj = fv | 0, a = fu - ii, b = fv - jj, k = jj * m.nu + ii, k2 = k + m.nu;
          px += (m.dx[k] * (1 - a) + m.dx[k + 1] * a) * (1 - b) + (m.dx[k2] * (1 - a) + m.dx[k2 + 1] * a) * b;
          py += (m.dy[k] * (1 - a) + m.dy[k + 1] * a) * (1 - b) + (m.dy[k2] * (1 - a) + m.dy[k2 + 1] * a) * b;
        }
        ax[i] = (px + 0.5) / sc - 0.5; ay[i] = (py + 0.5) / scY - 0.5;  // analysis → full-res pixel centres
      }
    }
  }

  return { detect, makeModel, fillMaps, stats, debug };
})();
