/**
 * OKLab / OKLCH color gradient splines (Unified)
 *
 * - method: 'natural-cubic' | 'centripetal-CR' | 'chordal-CR'  (controls tangent computation)
 * - param:  'uniform' | 'chordal' | 'centripetal'             (controls parameterization)
 *
 * Default: method='natural-cubic', param='uniform'
 */

/* ===================== Math Helpers ===================== */

const PI = Math.PI

const lin = (v) => (v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4)
const gam = (v) => (v <= 0.0031308 ? 12.92 * v : 1.055 * v ** (1 / 2.4) - 0.055)

const vScale = (v, s) => v.map((x) => x * s)
const vAdd = (a, b) => a.map((x, i) => x + (b[i] || 0))
const vSub = (a, b) => a.map((x, i) => x - (b[i] || 0))

// 3D distance excluding alpha channel (spatial geometry only)
function vDist3(a, b = []) {
  let sum = 0
  for (let i = 0; i < 3; i++) {
    const d = (a[i] || 0) - (b[i] || 0)
    sum += d * d
  }
  return Math.sqrt(sum)
}

/* ----------------- Color Conversions ----------------- */

function rgbToOklab(vec) {
  const [r0, g0, b0, alpha] = vec
  const r = lin(r0), g = lin(g0), b = lin(b0)
  const l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
  const m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
  const s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
  const l_ = Math.cbrt(l), m_ = Math.cbrt(m), s_ = Math.cbrt(s)
  const res = [
    0.2104542553 * l_ + 0.793617785 * m_ - 0.0040720468 * s_,
    1.9779984951 * l_ - 2.428592205 * m_ + 0.4505937099 * s_,
    0.0259040371 * l_ + 0.7827717662 * m_ - 0.808675766 * s_,
  ]
  if (alpha !== undefined) res.push(alpha)
  return res
}

function oklabToRgb(vec) {
  const [L, a, b, alpha] = vec
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b
  const s_ = L - 0.0894841775 * a - 1.291485548 * b
  const l = l_ ** 3, m = m_ ** 3, s = s_ ** 3
  const res = [
    gam(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s),
    gam(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s),
    gam(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s),
  ]
  if (alpha !== undefined) res.push(alpha)
  return res
}

function oklabToOklch(vec) {
  const [l, a, b, alpha] = vec
  const c = Math.hypot(a || 0, b || 0)
  const h = (Math.atan2(b || 0, a || 0) * 180) / PI
  const res = [l, c, (h + 360) % 360]
  if (alpha !== undefined) res.push(alpha)
  return res
}

function oklchToOklab(vec) {
  const [l, c, h, alpha] = vec
  const rad = (h * PI) / 180
  const res = [l, c * Math.cos(rad), c * Math.sin(rad)]
  if (alpha !== undefined) res.push(alpha)
  return res
}

/* ----------------- Parsing ----------------- */

const parseNumbers = (s) => (s && s.match(/[\d.+-eE]+/g) ? s.match(/[\d.+-eE]+/g).map(Number) : [])

const hexToRgb = (hex) => {
  let h = hex.replace('#', '').trim()
  if (h.length === 3) {
    const r = parseInt(h[0] + h[0], 16)
    const g = parseInt(h[1] + h[1], 16)
    const b = parseInt(h[2] + h[2], 16)
    return [r / 255, g / 255, b / 255]
  } else if (h.length === 4) {
    const r = parseInt(h[0] + h[0], 16)
    const g = parseInt(h[1] + h[1], 16)
    const b = parseInt(h[2] + h[2], 16)
    const a = parseInt(h[3] + h[3], 16)
    return [r / 255, g / 255, b / 255, a / 255]
  } else if (h.length === 6) {
    const n = parseInt(h, 16)
    return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255]
  } else if (h.length === 8) {
    const n = parseInt(h, 16)
    return [((n >> 24) & 255) / 255, ((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255]
  } else {
    const n = parseInt(h.slice(0, 6), 16)
    return [((n >> 16) & 255) / 255, ((n >> 8) & 255) / 255, (n & 255) / 255]
  }
}

const parseHsl = (s) => {
  const [h, sp, l, a] = parseNumbers(s)
  const res = [h, sp / 100, l / 100]
  if (a !== undefined) res.push(a)
  return res
}

const hslToRgb = (vec) => {
  const [h, s, l, a] = vec
  const f = (n) => {
    const k = (n + h / 30) % 12
    return l - s * Math.min(l, 1 - l) * Math.max(-1, Math.min(k - 3, 9 - k, 1))
  }
  const res = [f(0), f(8), f(4)]
  if (a !== undefined) res.push(a)
  return res
}

function parseColor(input, targetSpace = 'oklab') {
  if (Array.isArray(input)) return input.slice()
  if (typeof input !== 'string') throw new Error('Unsupported color input type')

  let lab
  if (input.startsWith('#')) lab = rgbToOklab(hexToRgb(input))
  else if (input.startsWith('rgb')) lab = rgbToOklab(parseNumbers(input).map((v, i) => (i < 3 ? v / 255 : v)))
  else if (input.startsWith('hsl')) lab = rgbToOklab(hslToRgb(parseHsl(input)))
  else if (input.startsWith('oklch')) {
    const lch = parseNumbers(input)
    if (targetSpace === 'oklch') return lch
    return oklchToOklab(lch)
  } else if (input.startsWith('oklab')) {
    const nums = parseNumbers(input)
    if (targetSpace === 'oklab') return nums
    return oklabToOklch(nums)
  } else throw new Error('Unsupported color: ' + input)

  return targetSpace === 'oklch' ? oklabToOklch(lab) : lab
}

/* ----------------- Hue Wrapping ----------------- */

function unwrapAngles(arr) {
  if (!arr || arr.length === 0) return []
  const out = arr.slice()
  for (let i = 1; i < out.length; i++) {
    let delta = out[i] - out[i - 1]
    while (delta <= -180) delta += 360
    while (delta > 180) delta -= 360
    out[i] = out[i - 1] + delta
  }
  return out
}

function wrapAngle(a) { return ((a % 360) + 360) % 360 }

function getDelta(pB, pA, space) {
  const d = pB.map((v, i) => v - (pA[i] || 0))
  if (space === 'oklch') {
    let dh = d[2]
    while (dh > 180) dh -= 360
    while (dh < -180) dh += 360
    d[2] = dh
  }
  return d
}

/* ===================== Solvers ===================== */

function solveTridiagonal(n, a, b, c, d) {
  const eps = 1e-14
  const cPrime = new Array(n)
  const dPrime = new Array(n)

  if (Math.abs(b[0]) < eps) b[0] = (b[0] >= 0 ? eps : -eps)
  cPrime[0] = c[0] / b[0]
  dPrime[0] = d[0] / b[0]

  for (let i = 1; i < n; i++) {
    let temp = b[i] - a[i] * cPrime[i - 1]
    if (Math.abs(temp) < eps) temp = (temp >= 0 ? eps : -eps)
    cPrime[i] = c[i] / temp
    dPrime[i] = (d[i] - a[i] * dPrime[i - 1]) / temp
  }

  const x = new Array(n)
  x[n - 1] = dPrime[n - 1]
  for (let i = n - 2; i >= 0; i--) {
    x[i] = dPrime[i] - cPrime[i] * x[i + 1]
  }
  return x
}

function solveCyclicTridiagonal(n, a, b, c, rhs) {
  if (n <= 2) return solveTridiagonal(n, a, b, c, rhs)

  let gamma = -b[0]
  if (Math.abs(gamma) < 1e-12) gamma = -1e-6

  const bb = [...b]
  bb[0] -= gamma
  bb[n - 1] -= (c[n - 1] * a[0]) / gamma

  const y = solveTridiagonal(n, a, bb, c, rhs)

  const u = new Array(n).fill(0)
  u[0] = gamma
  u[n - 1] = c[n - 1]
  const z = solveTridiagonal(n, a, bb, c, u)

  const denom = 1 + z[0] + (a[0] / gamma) * z[n - 1]
  if (Math.abs(denom) < 1e-14) {
    // Fallback for degenerate case
    for (let i = 0; i < bb.length; i++) bb[i] += 1e-12
    return solveTridiagonal(n, a, bb, c, rhs)
  }

  const fact = (y[0] + (a[0] / gamma) * y[n - 1]) / denom
  return y.map((val, i) => val - fact * z[i])
}

/* ===================== Spline Fitters ===================== */

function fitNaturalCubic(P, chords, loop, space) {
  const n = P.length
  const dims = Math.max(...P.map(p => p.length))
  const M = new Array(n).fill(0).map(() => new Array(dims).fill(0))

  for (let k = 0; k < dims; k++) {
    let vals = P.map(p => (p[k] !== undefined ? p[k] : 0))
    if (space === 'oklch' && k === 2) vals = unwrapAngles(vals)

    const a = new Array(n).fill(0)
    const b = new Array(n).fill(0)
    const c = new Array(n).fill(0)
    const r = new Array(n).fill(0)

    if (loop) {
      for (let i = 0; i < n; i++) {
        const prev = (i - 1 + n) % n
        const next = (i + 1) % n

        const hPrev = chords[prev]
        const hNext = chords[i]

        a[i] = hNext
        b[i] = 2 * (hNext + hPrev)
        c[i] = hPrev

        let valPrev = vals[prev]
        let valCurr = vals[i]
        let valNext = vals[next]

        // Handle hue wrapping at boundaries
        if (space === 'oklch' && k === 2) {
          if (i === 0) {
            let d = valPrev - valCurr
            while (d > 180) d -= 360
            while (d < -180) d += 360
            valPrev = valCurr + d
          }
          if (i === n - 1) {
            let d = valNext - valCurr
            while (d > 180) d -= 360
            while (d < -180) d += 360
            valNext = valCurr + d
          }
        }

        const slopePrev = (valCurr - valPrev) / hPrev
        const slopeNext = (valNext - valCurr) / hNext
        r[i] = 3 * (hNext * slopePrev + hPrev * slopeNext)
      }

      const res = solveCyclicTridiagonal(n, a, b, c, r)
      for (let i = 0; i < n; i++) M[i][k] = res[i]
    } else {
      for (let i = 1; i < n - 1; i++) {
        const hPrev = chords[i - 1]
        const hNext = chords[i]
        a[i] = hNext
        b[i] = 2 * (hNext + hPrev)
        c[i] = hPrev

        const slopePrev = (vals[i] - vals[i - 1]) / hPrev
        const slopeNext = (vals[i + 1] - vals[i]) / hNext
        r[i] = 3 * (hNext * slopePrev + hPrev * slopeNext)
      }

      // Natural boundary conditions
      const h0 = chords[0]
      b[0] = 2 * h0
      c[0] = h0
      r[0] = 3 * (vals[1] - vals[0])

      const hLast = chords[n - 2]
      a[n - 1] = hLast
      b[n - 1] = 2 * hLast
      r[n - 1] = 3 * (vals[n - 1] - vals[n - 2])

      const res = solveTridiagonal(n, a, b, c, r)
      for (let i = 0; i < n; i++) M[i][k] = res[i]
    }
  }

  return M
}

function fitCatmullRom(P, crChords, loop, space) {
  const n = P.length
  const M = new Array(n)

  for (let i = 0; i < n; i++) {
    // Endpoint handling for non-loop
    if (!loop && (i === 0 || i === n - 1)) {
      if (i === 0) M[i] = getDelta(P[1], P[0], space)
      else M[i] = getDelta(P[n - 1], P[n - 2], space)
      continue
    }

    const prevIdx = (i - 1 + n) % n
    const nextIdx = (i + 1) % n

    const segments = loop ? n : n - 1
    const dt0 = crChords[prevIdx % segments]  // chord from prev to i
    const dt1 = crChords[i % segments]         // chord from i to next

    const v1 = getDelta(P[i], P[prevIdx], space)
    const v2 = getDelta(P[nextIdx], P[i], space)

    const w1 = dt1 * dt1
    const w2 = dt0 * dt0
    const denom = Math.max(1e-12, dt0 * dt1 * (dt0 + dt1))

    M[i] = v1.map((val, k) => ((w1 * val + w2 * v2[k]) / denom) * dt1)
  }
  return M
}

/* ===================== Utilities ===================== */

function findSegmentIndex(u, t) {
  let lo = 0
  let hi = u.length - 2
  if (t <= u[0]) return 0
  if (t >= u[u.length - 1]) return u.length - 2
  while (lo <= hi) {
    const mid = Math.floor((lo + hi) / 2)
    if (t < u[mid]) hi = mid - 1
    else if (t >= u[mid + 1]) lo = mid + 1
    else return mid
  }
  return Math.max(0, Math.min(u.length - 2, lo))
}

/* ----------------- formatColor ----------------- */

export function formatColor(vec, format = 'oklab', inputSpace = 'oklab') {
  let lab = vec
  if (inputSpace === 'oklch') lab = oklchToOklab(vec)

  if (format === 'oklab') return lab
  if (format === 'oklch') return oklabToOklch(lab)

  const rgb = oklabToRgb(lab)
  const [r, g, b, a] = rgb

  if (format === 'rgb') {
    if (a !== undefined) return `rgba(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)},${+a.toFixed(3)})`
    return `rgb(${Math.round(r * 255)},${Math.round(g * 255)},${Math.round(b * 255)})`
  }

  // Hex format
  const to255 = (v) => Math.round(255 * Math.max(0, Math.min(1, v)))
  const hr = to255(r).toString(16).padStart(2, '0')
  const hg = to255(g).toString(16).padStart(2, '0')
  const hb = to255(b).toString(16).padStart(2, '0')
  if (a !== undefined) {
    const ha = to255(a).toString(16).padStart(2, '0')
    return `#${hr}${hg}${hb}${ha}`
  }
  return `#${hr}${hg}${hb}`
}

/* ===================== Main Logic ===================== */

export function fitSpline(fixpoints, { space = 'oklab', method = 'natural-cubic', param = 'uniform', loop = false, strengths = 1 } = {}) {
  if (!['oklab', 'oklch'].includes(space)) throw new Error('space must be "oklab" or "oklch"')
  if (!['natural-cubic', 'centripetal-CR', 'chordal-CR'].includes(method)) throw new Error('invalid method')
  if (!['uniform', 'chordal', 'centripetal'].includes(param)) throw new Error('invalid param option')

  const P = fixpoints.map((c) => parseColor(c, space))
  const n = P.length
  if (n < 2) throw new Error('Need at least 2 points')

  const segments = loop ? n : n - 1

  // Geometric chords (spatial distance only, excluding alpha)
  const geomChords = []
  for (let i = 0; i < segments; i++) {
    const a = P[i]
    const b = P[(i + 1) % n]
    const d = getDelta(b, a, space)
    let dist = vDist3(d)
    if (dist < 1e-6) dist = 1e-6
    geomChords.push(dist)
  }

  // Parameter chords: determine how 0..1 maps to fixpoints
  let paramChords = []
  if (param === 'uniform') {
    for (let i = 0; i < segments; i++) paramChords.push(1)
  } else if (param === 'chordal') {
    paramChords = geomChords.slice()
  } else { // centripetal
    paramChords = geomChords.map(d => Math.sqrt(d))
  }

  // CR chords: for tangent computation (geometric^power)
  const p = method === 'centripetal-CR' ? 0.5 : 1.0
  const crChords = geomChords.map(d => Math.pow(d, p))

  // Compute tangents
  let M
  if (method === 'natural-cubic') {
    // Returns derivatives D_i per unit paramChord
    M = fitNaturalCubic(P, paramChords, loop, space)
  } else {
    // Returns tangents in cr-chord units
    M = fitCatmullRom(P, crChords, loop, space)
  }

  // Apply strength scaling
  const getS = (i) => (Array.isArray(strengths) ? (strengths[i] ?? 1) : strengths)
  M = M.map((tan, i) => vScale(tan, getS(i)))

  // Build Hermite segments (each parameterized over t ∈ [0,1])
  const segs = []
  for (let i = 0; i < segments; i++) {
    const p0 = P[i]
    const p1 = P[(i + 1) % n]
    let m0 = M[i]
    let m1 = M[(i + 1) % n]

    // Convert tangents to Hermite form for this segment
    if (method === 'natural-cubic') {
      // Natural cubic returns derivatives; convert to Hermite tangents
      const h = paramChords[i]
      m0 = vScale(m0, h)
      m1 = vScale(m1, h)
    } else {
      // CR tangents are in cr-units; convert to param-units for this segment
      const safe = (x) => (Math.abs(x) < 1e-12 ? 1e-12 : x)
      m0 = vScale(m0, paramChords[i] / safe(crChords[i]))
      m1 = vScale(m1, paramChords[i] / safe(crChords[(i + 1) % segments]))
    }

    const dims = Math.max(p0.length, p1.length)
    const a = [], b = [], c = [], d = []

    // Handle OKLCH hue wrapping
    let targetP1 = p1
    if (space === 'oklch') {
      let dh = (p1[2] || 0) - (p0[2] || 0)
      while (dh > 180) dh -= 360
      while (dh < -180) dh += 360
      targetP1 = [...p1]
      targetP1[2] = (p0[2] || 0) + dh
    }

    // Build Hermite polynomial: P(t) = a + bt + ct² + dt³
    for (let k = 0; k < dims; k++) {
      const p0k = p0[k] !== undefined ? p0[k] : 0
      const p1k = targetP1[k] !== undefined ? targetP1[k] : 0
      const m0k = m0[k] !== undefined ? m0[k] : 0
      const m1k = m1[k] !== undefined ? m1[k] : 0
      a[k] = p0k
      b[k] = m0k
      const delta = p1k - p0k
      c[k] = 3 * delta - 2 * m0k - m1k
      d[k] = -2 * delta + m0k + m1k
    }

    segs.push({ a, b, c, d })
  }

  // Build cumulative parameter array
  const u = [0]
  for (let i = 0; i < segments; i++) u.push(u[i] + paramChords[i])

  return { segs, u, totalLen: u[u.length - 1], n, loop, space, method, param }
}

/**
 * Sample colors from a spline gradient
 * @param {Array|string[]} fixpoints - Color stops (arrays or CSS strings)
 * @param {number[]} weights - Sample positions in [0,1]
 * @param {Object} options - Configuration
 * @param {string} options.format - Output format: 'oklab'|'oklch'|'rgb'|'hex'
 * @param {string} options.space - Computation space: 'oklab'|'oklch'
 * @param {string} options.method - Spline method: 'natural-cubic'|'centripetal-CR'|'chordal-CR'
 * @param {string} options.param - Parameterization: 'uniform'|'chordal'|'centripetal'
 * @param {boolean} options.loop - Close the spline
 * @param {number|number[]} options.strengths - Tangent strength multipliers
 */
export function splineColors(fixpoints, weights, { format = 'oklab', space = 'oklab', method = 'natural-cubic', param = 'uniform', loop = false, strengths = 1 } = {}) {
  const spline = fitSpline(fixpoints, { space, method, param, loop, strengths })
  const u = spline.u

  return weights.map(w => {
    const clamped = Math.max(0, Math.min(1, w))
    const tGlobal = clamped * spline.totalLen

    const idx = findSegmentIndex(u, tGlobal)
    let i = idx
    if (tGlobal >= spline.totalLen) i = spline.segs.length - 1

    const len = u[i + 1] - u[i]
    const tLocal = len > 1e-9 ? (tGlobal - u[i]) / len : 0

    const s = spline.segs[i]
    const t2 = tLocal * tLocal
    const t3 = t2 * tLocal
    const res = []
    for (let k = 0; k < s.a.length; k++) {
      res[k] = s.a[k] + s.b[k] * tLocal + s.c[k] * t2 + s.d[k] * t3
    }

    if (space === 'oklch' && res.length >= 3) res[2] = wrapAngle(res[2])

    return formatColor(res, format, space)
  })
}