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
  function cholSolve(A, n, rhs) {          // A: dense SPD n×n (Float64Array, row-major), overwritten
    for (let j = 0; j < n; j++) {
      let d = A[j * n + j];
      for (let k = 0; k < j; k++) d -= A[j * n + k] ** 2;
      if (d <= 0) return null;
      d = Math.sqrt(d); A[j * n + j] = d;
      for (let i = j + 1; i < n; i++) {
        let s = A[i * n + j];
        for (let k = 0; k < j; k++) s -= A[i * n + k] * A[j * n + k];
        A[i * n + j] = s / d;
      }
    }
    return rhs.map(b => {
      const y = Float64Array.from(b);
      for (let i = 0; i < n; i++) { for (let k = 0; k < i; k++) y[i] -= A[i * n + k] * y[k]; y[i] /= A[i * n + i]; }
      for (let i = n - 1; i >= 0; i--) { for (let k = i + 1; k < n; k++) y[i] -= A[k * n + i] * y[k]; y[i] /= A[i * n + i]; }
      return y;
    });
  }
  function fitMesh(base, pts) {
    let u0 = Infinity, u1 = -Infinity, v0 = Infinity, v1 = -Infinity;
    for (const m of pts) { u0 = Math.min(u0, m.u); u1 = Math.max(u1, m.u); v0 = Math.min(v0, m.v); v1 = Math.max(v1, m.v); }
    const step = Math.max(4, Math.max(u1 - u0, v1 - v0) / 24);
    u0 -= step; v0 -= step;
    const nu = Math.ceil((u1 - u0) / step) + 2, nv = Math.ceil((v1 - v0) / step) + 2, n = nu * nv;
    const A = new Float64Array(n * n), bx = new Float64Array(n), by = new Float64Array(n);
    const mesh = { u0, v0, step, nu, nv };
    for (const m of pts) {
      const p = base.fwd0(m.u, m.v), rx = m.x - p[0], ry = m.y - p[1];
      const fu = (m.u - u0) / step, fv = (m.v - v0) / step, i = Math.floor(fu), j = Math.floor(fv), a = fu - i, b = fv - j;
      const ks = [j * nu + i, j * nu + i + 1, (j + 1) * nu + i, (j + 1) * nu + i + 1];
      const ws = [(1 - a) * (1 - b), a * (1 - b), (1 - a) * b, a * b];
      for (let x = 0; x < 4; x++) {
        bx[ks[x]] += ws[x] * rx; by[ks[x]] += ws[x] * ry;
        for (let y = 0; y < 4; y++) A[ks[x] * n + ks[y]] += ws[x] * ws[y];
      }
    }
    const alpha = 3, beta = 1e-3;          // membrane smoothness between neighbouring nodes, weak pull to zero
    for (let j = 0; j < nv; j++) for (let i = 0; i < nu; i++) {
      const k = j * nu + i; A[k * n + k] += beta;
      for (const kk of [i + 1 < nu ? k + 1 : -1, j + 1 < nv ? k + nu : -1]) {
        if (kk < 0) continue;
        A[k * n + k] += alpha; A[kk * n + kk] += alpha; A[k * n + kk] -= alpha; A[kk * n + k] -= alpha;
      }
    }
    const sol = cholSolve(A, n, [bx, by]);
    if (!sol) return null;
    mesh.dx = sol[0]; mesh.dy = sol[1];
    return mesh;
  }
  const rmsOf = (model, pts) => {
    let se = 0; for (const m of pts) { const p = model.forward(m.u, m.v); se += (p[0] - m.x) ** 2 + (p[1] - m.y) ** 2; }
    return Math.sqrt(se / pts.length);
  };
  // add the mesh only if it predicts held-out crossings better (2-fold check)
  function withMesh(model, inl) {
    if (inl.length < 80) return model;
    const A = inl.filter((m, k) => k % 2 === 0), B = inl.filter((m, k) => k % 2 === 1);
    let before = 0, after = 0;
    for (const [fitSet, testSet] of [[A, B], [B, A]]) {
      const mesh = fitMesh(model, fitSet);
      if (!mesh) return model;
      const mm = makeModel(model.H, model.lam, model.c, model.R, mesh);
      before += rmsOf(model, testSet); after += rmsOf(mm, testSet);
    }
    if (after > 0.95 * before) return model;
    const mesh = fitMesh(model, inl);
    const out = makeModel(model.H, model.lam, model.c, model.R, mesh);
    out.rmsPx = rmsOf(out, inl); out.inliers = inl; out.meshGain = before / after;
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
  function refineCrossing(t, W, H, q, e1, e2) {
    const L1 = Math.hypot(...e1), L2 = Math.hypot(...e2);
    const d1 = [e1[0] / L1, e1[1] / L1], d2 = [e2[0] / L2, e2[1] / L2];
    const n1 = [-d1[1], d1[0]], n2 = [-d2[1], d2[0]];
    const fit = (d, n, along, across) => {
      const S = [], O = [], h = 0.3 * across;
      for (const f of [-0.45, -0.38, -0.31, 0.31, 0.38, 0.45]) {
        const sa = f * along, bx = q[0] + d[0] * sa, by = q[1] + d[1] * sa;
        const vals = [];
        for (let o = -h; o <= h + 1e-9; o += 0.5) vals.push([o, bil(t, W, H, bx + n[0] * o, by + n[1] * o)]);
        let mn = Infinity; for (const v of vals) mn = Math.min(mn, v[1]);
        let sw = 0, so = 0;
        for (const [o, v] of vals) { const w = (v - mn) ** 2; sw += w; so += w * o; }
        if (sw > 0) { S.push(sa); O.push(so / sw); }
      }
      if (S.length < 4) return null;
      let ms = 0, mo = 0; for (let k = 0; k < S.length; k++) { ms += S[k]; mo += O[k]; }
      ms /= S.length; mo /= S.length;
      let sso = 0, sss = 0; for (let k = 0; k < S.length; k++) { sso += (S[k] - ms) * (O[k] - mo); sss += (S[k] - ms) ** 2; }
      const bb = sss ? sso / sss : 0;
      return { a: mo - bb * ms, b: bb };
    };
    const A = fit(d1, n1, L1, L2), B = fit(d2, n2, L2, L1);
    if (!A || !B) return null;
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
  function buildLattice(pts, W, H, centre) {
    const cell = 24, cols = Math.ceil(W / cell), rows = Math.ceil(H / cell);
    const buckets = new Map();
    pts.forEach((p, i) => {
      const k = Math.floor(p[1] / cell) * cols + Math.floor(p[0] / cell);
      if (!buckets.has(k)) buckets.set(k, []); buckets.get(k).push(i);
    });
    const used = new Uint8Array(pts.length);
    function nearest(x, y, r, skipUsed = true) {
      let best = -1, bd = r;
      const c0 = Math.floor((x - r) / cell), c1 = Math.floor((x + r) / cell), r0 = Math.floor((y - r) / cell), r1 = Math.floor((y + r) / cell);
      for (let rr = Math.max(0, r0); rr <= Math.min(rows - 1, r1); rr++)
        for (let cc = Math.max(0, c0); cc <= Math.min(cols - 1, c1); cc++) {
          const b = buckets.get(rr * cols + cc); if (!b) continue;
          for (const i of b) {
            if (skipUsed && used[i]) continue;
            const d = Math.hypot(pts[i][0] - x, pts[i][1] - y);
            if (d < bd) { bd = d; best = i; }
          }
        }
      return best;
    }
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
    const has = (x, y, tol) => nearest(x, y, tol, false) >= 0;

    // seed: a point whose 3×3 neighbourhood forms a clean square lattice
    const order = pts.map((p, i) => i).sort((a, b) => dist(pts[a], centre) - dist(pts[b], centre));
    let seed = -1, e1, e2;
    for (const i of order.slice(0, 1500)) {
      const p = pts[i], nb = neighbours(p, 5, 90);
      if (nb.length < 4) continue;
      const q1 = pts[nb[0][0]], a = [q1[0] - p[0], q1[1] - p[1]], la = Math.hypot(...a);
      let b = null;
      for (const [j, d] of nb) {
        if (d > 1.5 * la || d < 0.66 * la) continue;
        const v = [pts[j][0] - p[0], pts[j][1] - p[1]];
        if (Math.abs((v[0] * a[0] + v[1] * a[1]) / (d * la)) < 0.35) { b = v; break; }
      }
      if (!b) continue;
      const tol = 0.22 * la;
      const ok = has(p[0] - a[0], p[1] - a[1], tol) && has(p[0] - b[0], p[1] - b[1], tol) &&
                 has(p[0] + a[0] + b[0], p[1] + a[1] + b[1], tol) && has(p[0] - a[0] - b[0], p[1] - a[1] - b[1], tol) &&
                 has(p[0] + a[0] - b[0], p[1] + a[1] - b[1], tol) && has(p[0] - a[0] + b[0], p[1] - a[1] + b[1], tol) &&
                 has(p[0] + 2 * a[0], p[1] + 2 * a[1], tol * 1.5) && has(p[0] + 2 * b[0], p[1] + 2 * b[1], tol * 1.5);
      if (!ok) continue;
      seed = i; e1 = a; e2 = b; break;
    }
    if (seed < 0) return null;
    // orient: e1 ≈ image +x, e2 ≈ image +y
    if (Math.abs(e1[0]) < Math.abs(e2[0])) [e1, e2] = [e2, e1];
    if (e1[0] < 0) e1 = [-e1[0], -e1[1]];
    if (e2[1] < 0) e2 = [-e2[0], -e2[1]];

    const nodes = new Map(), key = (i, j) => i + ',' + j;
    const queue = [];
    const add = (i, j, pi, a, b) => { used[pi] = 1; const n = { i, j, x: pts[pi][0], y: pts[pi][1], e1: a, e2: b }; nodes.set(key(i, j), n); queue.push(n); };
    add(0, 0, seed, e1, e2);
    while (queue.length && nodes.size < 30000) {
      const n = queue.shift();
      for (const [di, dj] of [[1, 0], [-1, 0], [0, 1], [0, -1]]) {
        const ti = n.i + di, tj = n.j + dj;
        if (nodes.has(key(ti, tj))) continue;
        const back = nodes.get(key(n.i - di, n.j - dj));
        const base = di ? [n.e1[0] * di, n.e1[1] * di] : [n.e2[0] * dj, n.e2[1] * dj];
        let e = back ? [n.x - back.x, n.y - back.y] : base;
        if (Math.abs(Math.hypot(...e) / Math.hypot(...base) - 1) > 0.3) e = base;
        const le = Math.hypot(...e);
        const k = nearest(n.x + e[0], n.y + e[1], 0.3 * le);
        if (k < 0) continue;
        const step = [(pts[k][0] - n.x), (pts[k][1] - n.y)];
        const ne1 = di ? [step[0] * di, step[1] * di] : n.e1, ne2 = dj ? [step[0] * dj, step[1] * dj] : n.e2;
        add(ti, tj, k, ne1, ne2);
      }
    }
    return [...nodes.values()].map(n => ({ u: n.i, v: n.j, x: n.x, y: n.y }));
  }

  // ---------- predicted-position sweep with the fitted model ----------
  function sweep(model, t, W, H, uvBox, thrT) {
    const out = [];
    for (let v = uvBox.v0; v <= uvBox.v1; v++) for (let u = uvBox.u0; u <= uvBox.u1; u++) {
      const p = model.forward(u, v);
      if (!(p[0] > 2 && p[1] > 2 && p[0] < W - 3 && p[1] < H - 3)) continue;
      const pu = model.forward(u + 1, v), pv = model.forward(u, v + 1);
      const e1 = [pu[0] - p[0], pu[1] - p[1]], e2 = [pv[0] - p[0], pv[1] - p[1]];
      const s = Math.min(Math.hypot(...e1), Math.hypot(...e2));
      if (!(s > 4)) continue;
      const q = centroid(t, W, H, p[0], p[1], Math.max(2, Math.round(0.3 * s)));
      if (!q || Math.hypot(q[0] - p[0], q[1] - p[1]) > 0.25 * s) continue;
      // a real crossing: bright along all 4 arms, dark in the middle of the 4 neighbouring squares
      const arm = (dx, dy) => sampleMax(t, W, H, q[0] + dx, q[1] + dy) >= thrT;
      if (!(arm(0.4 * e1[0], 0.4 * e1[1]) && arm(-0.4 * e1[0], -0.4 * e1[1]) && arm(0.4 * e2[0], 0.4 * e2[1]) && arm(-0.4 * e2[0], -0.4 * e2[1]))) continue;
      let dark = 0;
      for (const [a, b] of [[0.5, 0.5], [-0.5, 0.5], [0.5, -0.5], [-0.5, -0.5]]) {
        const x = Math.round(q[0] + a * e1[0] + b * e2[0]), y = Math.round(q[1] + a * e1[1] + b * e2[1]);
        if (x >= 0 && y >= 0 && x < W && y < H && t[y * W + x] < thrT) dark++;
      }
      if (dark < 3) continue;
      const fine = refineCrossing(t, W, H, q, e1, e2);
      if (!fine) continue;
      out.push({ u, v, x: fine[0], y: fine[1] });
    }
    return out;
  }

  // ---------- main entry ----------
  // bitmap → { model (analysis coords), sc, W, H, matches, rmsMm, spacingPx, visible:{u0,u1,v0,v1} }
  function detect(bitmap, maskFn, maxDim = 1600) {
    const scale = Math.min(1, maxDim / Math.max(bitmap.width, bitmap.height));
    const c = document.createElement('canvas');
    c.width = Math.round(bitmap.width * scale); c.height = Math.round(bitmap.height * scale);
    c.getContext('2d').drawImage(bitmap, 0, 0, c.width, c.height);
    const sc = c.width / bitmap.width, W = c.width, H = c.height;
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
    cv.findContours(mask, contours, hier, cv.RETR_EXTERNAL, cv.CHAIN_APPROX_SIMPLE);
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

    const fail = reason => ({ ok: false, reason, sc, W, H, candidates: pts.length });
    if (pts.length < 30 || !outline.length) return fail('no grid lines found on the green area');
    let mx = 0, my = 0; for (const p of outline) { mx += p[0]; my += p[1]; }
    const lattice = buildLattice(pts, W, H, [mx / outline.length, my / outline.length]);
    if (!lattice || lattice.length < 30) return fail('could not lock onto the 1 cm grid');
    // only green regions that actually carry the grid count as mat (ignores other green objects)
    const onMat = regions.filter(poly => lattice.filter(m => inPoly(poly, m.x, m.y)).length >= 5);
    if (onMat.length) outline = onMat.flat();

    const cImg = [W / 2, H / 2], R = Math.hypot(W, H) / 2;
    let model = fitModel(lattice, cImg, R, false);
    if (!model) return fail('grid fit failed');
    // thresholds for crossing validation: typical line brightness
    const vals = model.inliers.map(m => t[Math.round(m.y) * W + Math.round(m.x)]).sort((a, b) => a - b);
    const thrT = 0.35 * vals[vals.length >> 1];
    let matches = lattice;
    for (let round = 0; round < 3; round++) {
      const box = uvBounds(model, outline, matches, 2);
      const sw = sweep(model, t, W, H, box, thrT);
      if (sw.length < 30) break;
      matches = sw;
      const m2 = fitModel(matches, cImg, R, true);
      if (!m2) break;
      model = m2;
    }
    model = withMesh(model, model.inliers);
    const inl = model.inliers;
    // error in cm space
    let se = 0; const errs = [];
    for (const m of inl) { const q = model.inverse(m.x, m.y); const e = Math.hypot(q[0] - m.u, q[1] - m.v); se += e * e; errs.push(e); }
    errs.sort((a, b) => a - b);
    const visible = uvBounds(model, outline, inl, 0);
    return {
      ok: true, model, sc, W, H, matches: inl, candidates: pts.length,
      rmsMm: Math.sqrt(se / inl.length) * 10, p95Mm: errs[Math.floor(errs.length * 0.95)] * 10,
      spacingPx: medianSpacing(inl) / sc, visible, lam: model.lam, outline, mesh: !!model.corr,
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
  function uvBounds(model, outline, matches, pad) {
    let u0 = Infinity, u1 = -Infinity, v0 = Infinity, v1 = -Infinity;
    for (const p of outline) {
      const q = model.inverse(p[0], p[1]);
      if (!isFinite(q[0]) || !isFinite(q[1])) continue;
      u0 = Math.min(u0, q[0]); u1 = Math.max(u1, q[0]); v0 = Math.min(v0, q[1]); v1 = Math.max(v1, q[1]);
    }
    let mu0 = Infinity, mu1 = -Infinity, mv0 = Infinity, mv1 = -Infinity;
    for (const m of matches) { mu0 = Math.min(mu0, m.u); mu1 = Math.max(mu1, m.u); mv0 = Math.min(mv0, m.v); mv1 = Math.max(mv1, m.v); }
    const lim = 15;
    return {
      u0: Math.floor(Math.max(u0, mu0 - lim)) - pad, u1: Math.ceil(Math.min(u1, mu1 + lim)) + pad,
      v0: Math.floor(Math.max(v0, mv0 - lim)) - pad, v1: Math.ceil(Math.min(v1, mv1 + lim)) + pad,
    };
  }

  // output pixel grid → analysis-image coordinates, written into ax/ay (Float32Array) for rows y0..y0+rows
  function fillMaps(model, r, S, y0, rows, W, ax, ay, sc) {
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
        ax[i] = (px + 0.5) / sc - 0.5; ay[i] = (py + 0.5) / sc - 0.5;   // analysis → full-res pixel centres
      }
    }
  }

  return { detect, makeModel, fillMaps };
})();
