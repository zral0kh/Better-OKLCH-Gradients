// Refs: Catmull & Rom 1974; Yuksel et al. 2011 (Parameterization)
// now also with chordal catmull-rom for better color spacing

//////////////////// MATH HELPERS ////////////////////

const PI = Math.PI
const TAU = 2 * PI

// Linear <-> sRGB Gamma
const lin = (v) => (v <= 0.04045 ? v / 12.92 : ((v + 0.055) / 1.055) ** 2.4)
const gam = (v) => (v <= 0.0031308 ? 12.92 * v : 1.055 * v ** (1 / 2.4) - 0.055)

// Euclidean Norm (N-dimensional)
const vLen = (v) => Math.sqrt(v.reduce((sum, x) => sum + x * x, 0))

// Vector Scale
const vScale = (v, s) => v.map((x) => x * s)

// Calculate Delta vector between two points, handling Hue wrapping if in Polar space
// In OKLCH arrays: [L, C, H] -> index 2 is Hue
function getDelta(pB, pA, space) {
  const d = pB.map((v, i) => v - pA[i])
  if (space === 'oklch') {
    let dh = d[2]
    // Wrap hue to shortest path (-180 to 180)
    while (dh > 180) dh -= 360
    while (dh < -180) dh += 360
    d[2] = dh
  }
  return d
}

//////////////////// COLOR CONVERSION ////////////////////

// RGB -> OKLAB
function rgbToOklab([r, g, b]) {
  r = lin(r)
  g = lin(g)
  b = lin(b)
  const l = 0.4122214708 * r + 0.5363325363 * g + 0.0514459929 * b
  const m = 0.2119034982 * r + 0.6806995451 * g + 0.1073969566 * b
  const s = 0.0883024619 * r + 0.2817188376 * g + 0.6299787005 * b
  const l_ = Math.cbrt(l), m_ = Math.cbrt(m), s_ = Math.cbrt(s)
  return [
    0.2104542553 * l_ + 0.793617785 * m_ - 0.0040720468 * s_,
    1.9779984951 * l_ - 2.428592205 * m_ + 0.4505937099 * s_,
    0.0259040371 * l_ + 0.7827717662 * m_ - 0.808675766 * s_,
  ]
}

// OKLAB -> RGB
function oklabToRgb([L, a, b]) {
  const l_ = L + 0.3963377774 * a + 0.2158037573 * b
  const m_ = L - 0.1055613458 * a - 0.0638541728 * b
  const s_ = L - 0.0894841775 * a - 1.291485548 * b
  const l = l_ ** 3, m = m_ ** 3, s = s_ ** 3
  return [
    gam(4.0767416621 * l - 3.3077115913 * m + 0.2309699292 * s),
    gam(-1.2684380046 * l + 2.6097574011 * m - 0.3413193965 * s),
    gam(-0.0041960863 * l - 0.7034186147 * m + 1.707614701 * s),
  ]
}

// OKLAB -> OKLCH
function oklabToOklch([l, a, b]) {
  const c = Math.hypot(a, b)
  const h = (Math.atan2(b, a) * 180) / PI
  return [l, c, (h + 360) % 360]
}

// OKLCH -> OKLAB
function oklchToOklab([l, c, h]) {
  const rad = (h * PI) / 180
  return [l, c * Math.cos(rad), c * Math.sin(rad)]
}

// Generic Parser
function parseColor(c, targetSpace = 'oklab') {
  let lab
  if (Array.isArray(c)) lab = c // Assume input is already [L,a,b] or [L,c,h] if raw array
  else if (typeof c === 'string') {
    if (c.startsWith('#')) lab = rgbToOklab(hexToRgb(c))
    else if (c.startsWith('rgb')) lab = rgbToOklab(parseNumbers(c).map(v => v/255))
    else if (c.startsWith('hsl')) lab = rgbToOklab(hslToRgb(parseHsl(c)))
    else if (c.startsWith('oklch')) {
      const lch = parseNumbers(c) // [L, C, H]
      lab = targetSpace === 'oklch' ? lch : oklchToOklab(lch)
      if (targetSpace === 'oklch') return lab // Short circuit
    } else {
        throw new Error('Unsupported color format: ' + c)
    }
  }
  
  // Convert LAB result to LCH if requested
  if (targetSpace === 'oklch' && !c.startsWith('oklch')) {
    return oklabToOklch(lab)
  }
  return lab
}

// String Helpers
const parseNumbers = (s) => s.match(/[\d.-]+/g).map(Number)
const hexToRgb = (hex) => {
  const n = parseInt(hex.slice(1), 16)
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255].map((v) => v / 255)
}
const parseHsl = (s) => {
  const [h, s_, l] = parseNumbers(s)
  return [h, s_ / 100, l / 100]
}
const hslToRgb = ([h, s, l]) => {
  const f = (n) => {
    const k = (n + h / 30) % 12
    return l - s * Math.min(l, 1 - l) * Math.max(-1, Math.min(k - 3, 9 - k, 1))
  }
  return [f(0), f(8), f(4)]
}

//////////////////// SPLINE LOGIC ////////////////////

/**
 * Fits a Cubic Hermite Spline using Non-Uniform (Chordal) Catmull-Rom parameterization.
 * 
 * @param {string[]} fixpoints - Array of color strings.
 * @param {Object} options - Configuration object.
 * @param {string} [options.space='oklab'] - Computation space: 'oklab' (Cartesian) or 'oklch' (Polar).
 * @param {boolean} [options.loop=false] - Whether to close the loop (last point connects to first).
 * @param {number|number[]} [options.strengths=1] - Tension scalar(s).
 * @returns {Object} - Spline data { P, M, n, space, loop }.
 */
export function fitSpline(fixpoints, { space = 'oklab', loop = false, strengths = 1 } = {}) {
  const P = fixpoints.map(c => parseColor(c, space))
  const n = P.length
  if (n < 2) throw new Error("Need at least 2 points")

  // Resolve strength scalar/array
  const getS = (i) => (Array.isArray(strengths) ? (strengths[i] ?? 1) : strengths)

  // 1. Calculate Chord Lengths (Euclidean distance in selected space)
  const chords = []
  const segments = loop ? n : n - 1
  
  for (let i = 0; i < segments; i++) {
    const pA = P[i]
    const pB = P[(i + 1) % n]
    let dist = vLen(getDelta(pB, pA, space))
    if (dist < 0.0001) dist = 0.0001 // Avoid div by zero
    chords.push(dist)
  }

  // 2. Compute Tangents (M) using Chordal weights
  const M = new Array(n)

  for (let i = 0; i < n; i++) {
    let tangent

    if (!loop && (i === 0 || i === n - 1)) {
        // Open Curve Endpoints: Use Linear Projection (Lowest Energy for 2 points)
        if (i === 0) tangent = getDelta(P[1], P[0], space)
        else tangent = getDelta(P[n - 1], P[n - 2], space)
    } else {
        // Internal or Loop points: Non-Uniform Catmull-Rom
        const prevIdx = (i - 1 + n) % n
        const nextIdx = (i + 1) % n
        
        const dt0 = chords[(!loop && i===0) ? 0 : (prevIdx) % segments] // dist(prev, curr)
        const dt1 = chords[(!loop && i===n-1) ? segments-1 : i % segments] // dist(curr, next)

        const v1 = getDelta(P[i], P[prevIdx], space) // P_i - P_prev
        const v2 = getDelta(P[nextIdx], P[i], space) // P_next - P_i

        // Weighted Tangent Formula
        // T_norm = (dt1^2 * v1 + dt0^2 * v2) / (dt0 * dt1 * (dt0 + dt1))
        // We multiply by dt1 to normalize for the outgoing Hermite interval [0, 1]
        // M_i = T_norm * dt1
        const w1 = dt1 * dt1
        const w2 = dt0 * dt0
        const denom = dt0 * (dt0 + dt1) // dt1 canceled out from numerator scaling
        
        tangent = v1.map((val, k) => (w1 * val + w2 * v2[k]) / denom)
    }

    // Apply user strength/tension
    M[i] = vScale(tangent, getS(i))
  }

  return { P, M, n, space, loop, chords }
}

/**
 * Evaluates the spline at a specific segment and local t.
 */
export function evalHermite(p0, m0, p1, m1, t, space) {
  const t2 = t * t
  const t3 = t2 * t
  
  const h00 = 2 * t3 - 3 * t2 + 1
  const h10 = t3 - 2 * t2 + t
  const h01 = -2 * t3 + 3 * t2
  const h11 = t3 - t2

  // Handle Hue Wrapping for P1 relative to P0 if in Polar space
  let targetP1 = p1
  if (space === 'oklch') {
    let dh = p1[2] - p0[2]
    while (dh > 180) dh -= 360
    while (dh < -180) dh += 360
    if (Math.abs(dh - (p1[2] - p0[2])) > 0.001) {
        targetP1 = [p1[0], p1[1], p0[2] + dh]
    }
  }

  const result = p0.map((_, i) => 
    h00 * p0[i] + h10 * m0[i] + h01 * targetP1[i] + h11 * m1[i]
  )

  // Normalize Hue result
  if (space === 'oklch') {
    result[2] = (result[2] % 360 + 360) % 360
  }
  return result
}

//////////////////// OUTPUT & UTIL ////////////////////

export function formatColor(vec, format = 'oklab', inputSpace = 'oklab') {
  // Normalize input to OKLAB for conversion
  const lab = inputSpace === 'oklch' ? oklchToOklab(vec) : vec

  if (format === 'oklab') return lab
  if (format === 'oklch') {
    return inputSpace === 'oklch' 
        ? { l: vec[0], c: vec[1], h: vec[2] } 
        : { l: lab[0], c: Math.hypot(lab[1], lab[2]), h: (Math.atan2(lab[2], lab[1]) * 180 / PI + 360) % 360 }
  }
  
  const [r, g, b] = oklabToRgb(lab).map((v) => Math.round(255 * Math.min(1, Math.max(0, v))))
  
  if (format === 'hex') {
    return '#' + [r, g, b].map((v) => v.toString(16).padStart(2, '0')).join('')
  }
  return `rgb(${r},${g},${b})`
}

export function oklchToCss(oklch) {
  const L = Array.isArray(oklch) ? oklch[0] : oklch.l
  const C = Array.isArray(oklch) ? oklch[1] : oklch.c
  const H = Array.isArray(oklch) ? oklch[2] : oklch.h
  return `oklch(${Math.round(L*1000)/10}% ${Math.round(C*1000)/1000} ${Math.round(H*100)/100}deg)`
}

/**
 * Maps global 0..1 weights to segment indices and local t.
 * For Chordal splines, uniform parameterization along the array indices isn't perfectly time-accurate
 * relative to distance, but it allows predictable sampling count.
 */
export function prepareWeights(weights, n, loop) {
  const segments = loop ? n : n - 1
  return weights.map((w) => {
    // Clamp
    w = Math.max(0, Math.min(1, w))
    
    // Map to segment
    const scaled = w * segments
    let i = Math.floor(scaled)
    if (i >= segments) i = segments - 1
    
    return { i, t: scaled - i }
  })
}

/**
 * Generate interpolated colors.
 * @param {string[]} fixpoints - Array of CSS colors.
 * @param {number[]} weights - Array of numbers 0..1.
 * @param {Object} [options] - format, space, loop, strengths.
 */
export function splineColors(fixpoints, weights, { format = 'oklab', space = 'oklab', loop = false, strengths = 1 } = {}) {
  const spline = fitSpline(fixpoints, { space, loop, strengths })
  const prepped = prepareWeights(weights, spline.n, loop)
  
  return prepped.map(({ i, t }) => {
    // M[i] is already scaled for the segment STARTING at i.
    // However, if we are in Non-Uniform mode (Chordal), the 'incoming' tangent M[i+1] 
    // was scaled for the segment starting at i+1. We must rescale it for segment i.
    
    const p0 = spline.P[i]
    const p1 = spline.P[(i + 1) % spline.n]
    const m0 = spline.M[i]
    
    let m1 = spline.M[(i + 1) % spline.n]
    
    // Rescale M1 for the current segment length if using Chordal logic (implicit in fitSpline)
    if (spline.chords) {
        const distCurr = spline.chords[i]
        const distNext = spline.chords[(i + 1) % (loop ? spline.n : spline.n - 1)] || 1
        const scale = (distNext > 0.0001) ? (distCurr / distNext) : 1
        m1 = vScale(m1, scale)
    }

    const val = evalHermite(p0, m0, p1, m1, t, space)
    return formatColor(val, format, space)
  })
}