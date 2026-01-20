// Minimal spline-based OKLab interpolation (fit once, evaluate fast)
// Refs: Catmull & Rom 1974; Bartels et al. 1987; Ottosson 2020

//////////////////// COLOR PARSING ////////////////////

function parseColor(c) {
  if (typeof c === 'string') {
    if (c.startsWith('#')) return rgbToOklab(hexToRgb(c))
    if (c.startsWith('rgb')) return rgbToOklab(parseRgb(c))
    if (c.startsWith('hsl')) return rgbToOklab(hslToRgb(parseHsl(c)))
    if (c.startsWith('oklch')) return parseOklch(c)
  }
  throw new Error('Unsupported color format')
}

function parseOklch(s) {
  const [l, c, h] = s.match(/[\d.]+/g).map(Number)
  const rad = (h * Math.PI) / 180
  return [l, c * Math.cos(rad), c * Math.sin(rad)]
}

function hexToRgb(hex) {
  const n = parseInt(hex.slice(1), 16)
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255].map((v) => v / 255)
}

function parseRgb(s) {
  return s.match(/[\d.]+/g).map((v) => +v / 255)
}

function parseHsl(s) {
  const [h, s_, l] = s.match(/[\d.]+/g).map(Number)
  return [h, s_ / 100, l / 100]
}

function hslToRgb([h, s, l]) {
  const f = (n) => {
    const k = (n + h / 30) % 12
    return l - s * Math.min(l, 1 - l) * Math.max(-1, Math.min(k - 3, 9 - k, 1))
  }
  return [f(0), f(8), f(4)]
}

//////////////////// OKLAB ////////////////////

const lin = (v) => (v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4)
const gam = (v) => (v <= 0.0031308 ? 12.92 * v : 1.055 * v ** (1 / 2.4) - 0.055)

function rgbToOklab([r, g, b]) {
  r = lin(r)
  g = lin(g)
  b = lin(b)
  const l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
  const m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
  const s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
  const l_ = Math.cbrt(l),
    m_ = Math.cbrt(m),
    s_ = Math.cbrt(s)
  return [
    0.2104542553 * l_ + 0.793617785 * m_ - 0.0040720468 * s_,
    1.9779984951 * l_ - 2.428592205 * m_ + 0.4505937099 * s_,
    0.0259040371 * l_ + 0.7827717662 * m_ - 0.808675766 * s_,
  ]
}

function oklabToRgb([L, a, b]) {
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b
  const s_ = L - 0.0894841775 * a - 1.291485548 * b
  const l = l_ ** 3,
    m = m_ ** 3,
    s = s_ ** 3
  return [
    gam(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s),
    gam(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s),
    gam(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s),
  ]
}

//////////////////// SPLINE (FIT ONCE) ////////////////////

function computeTangents(P, strengths = 1) {
  const n = P.length
  // Resolve strength: if array, pick index; if number, use global; else 1
  const getS = (i) => (Array.isArray(strengths) ? (strengths[i] ?? 1) : strengths)

  return P.map((_, i) => {
    const s = getS(i)
    let vec
    
    // Standard Catmull-Rom logic (0.5 * finite difference)
    // We apply the 'strength' scalar to the result.
    if (i === 0) {
      // Start: Forward difference
      vec = P[1].map((v, j) => v - P[0][j])
    } else if (i === n - 1) {
      // End: Backward difference
      vec = P[n - 1].map((v, j) => v - P[n - 2][j])
    } else {
      // Middle: Central difference * 0.5
      vec = P[i + 1].map((v, j) => 0.5 * (v - P[i - 1][j]))
    }

    return vec.map((v) => v * s)
  })
}

/**
 * Fits a Cubic Hermite Spline (Catmull-Rom) to the provided colors.
 * @param {string[]} fixpoints - Array of CSS color strings.
 * @param {number|number[]} [strengths=1] - Tangent strength(s). 0 = sharp/linear, 1 = smooth.
 * @returns {Object} - Spline representation with control points P, tangents M, and number of points n.
 */
export function fitSpline(fixpoints, strengths = 1) {
  const P = fixpoints.map(parseColor)
  const M = computeTangents(P, strengths)
  return { P, M, n: P.length }
}

/**
 * Evaluates the Cubic Hermite Spline at parameter t in [0, 1].
 * @param {Array} p0 - Start point.
 * @param {Array} m0 - Start tangent.
 * @param {Array} p1 - End point.
 * @param {Array} m1 - End tangent.
 * @param {number} t - Parameter in [0, 1].
 * @returns {Array} - Interpolated point.
 */
export function evalHermite(p0, m0, p1, m1, t) {
  const t2 = t * t,
    t3 = t2 * t
  return p0.map(
    (_, i) =>
      (2 * t3 - 3 * t2 + 1) * p0[i] + (t3 - 2 * t2 + t) * m0[i] + (-2 * t3 + 3 * t2) * p1[i] + (t3 - t2) * m1[i],
  )
}

/**
 * 
 * @param {Array} lab - oklab color as [L, a, b].
 * @param {string} format - Output format: 'oklab' (array), 'oklch' (object), 'rgb' (CSS rgb string), 'hex' (hex string).
 * @returns {Array|string|Object} - Formatted color.
 */
export function formatColor(lab, format = 'oklab') {
  if (format === 'oklab') return lab
  if (format === 'oklch') {
    const [L, a, b] = lab
    return {
      l: L,
      c: Math.hypot(a, b),
      h: ((Math.atan2(b, a) * 180) / Math.PI + 360) % 360,
    }
  }
  const [r, g, b] = oklabToRgb(lab).map((v) => Math.round(255 * Math.min(1, Math.max(0, v))))
  if (format === 'hex') return '#' + [r, g, b].map((v) => v.toString(16).padStart(2, '0')).join('')
  return `rgb(${r},${g},${b})`
}

/**
 * 
 * @param {Array|Object} oklch - OKLCH color as array [L, C, H] or object {l, c, h}.
 * @returns {string} - CSS OKLCH color string.
 */
export function oklchToCss(oklch) {
  const L = Array.isArray(oklch) ? oklch[0] : oklch.l
  const C = Array.isArray(oklch) ? oklch[1] : oklch.c
  const H = Array.isArray(oklch) ? oklch[2] : oklch.h
  const lStr = `${Math.round(L * 1000) / 10}%`
  const cStr = `${Math.round(C * 1000) / 1000}`
  const hStr = `${Math.round(H * 100) / 100}deg`
  return `oklch(${lStr} ${cStr} ${hStr})`
}

/**
 * 
 * @param {number[]} weights - Array of numbers 0..1 to sample.
 * @param {number} n - Number of segments in the spline. If not provided, chooses based on weights.
 * @returns {Array} - Array of objects {i, t} where i is segment index and t is local parameter.
 */

export function prepareWeights(weights, n=-1) {
  if (n==-1) n = weights.length
  const maxSeg = n - 1
  return weights
    .map((w) => Math.min(1, Math.max(0, w)))
    .sort((a, b) => a - b)
    .map((w) => {
      const x = w * maxSeg
      const i = Math.min(maxSeg - 1, Math.floor(x))
      return { i, t: x - i }
    })
}

/**
 * Generate interpolated colors from a set of fixpoints.
 * @param {string[]} fixpoints - Array of CSS colors.
 * @param {number[]} weights - Array of numbers 0..1 to sample.
 * @param {string} [format='oklab'] - Output format ('oklab', 'oklch', 'rgb', 'hex').
 * @param {number|number[]} [strengths=1] - Tangent strength(s).
 * @returns {Array} - Array of interpolated colors in the specified format.
 */
export function splineColors(fixpoints, weights, format = 'oklab', strengths = 1) {
  const spline = fitSpline(fixpoints, strengths) // O(n) once
  const preppedWeights = prepareWeights(weights, spline.n) // O(m) once
  return preppedWeights.map((w) =>
    formatColor(evalHermite(spline.P[w.i], spline.M[w.i], spline.P[w.i + 1], spline.M[w.i + 1], w.t), format),
  )
}