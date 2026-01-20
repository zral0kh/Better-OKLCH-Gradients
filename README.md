# Better OKLCH Spline Gradient

This project provides a standalone JavaScript implementation and a web-based visualization tool for generating **perceptually smooth, high-order color gradients**.

By utilizing **Cubic Hermite Splines (Catmull-Rom)** within the **OKLCH** color space, this approach mitigates common artifacts found in standard linear sRGB interpolation, such as Mach bands and desaturated intermediate colors.

[![Live Demo Preview](assets/preview.png)](https://zral0kh.github.io/hermite-oklch-gpage/)

[**🚀 Open Live Editor**](https://zral0kh.github.io/hermite-oklch-gpage/)

---

## 1. Limitations of Standard Linear Gradients

Standard CSS gradients (`linear-gradient`) typically rely on **Linear Interpolation (Lerp)** within the **sRGB** color space. While computationally efficient, this method introduces three significant perceptual flaws:

### A. The "Gray Dead Zone"
When interpolating between complementary colors (e.g., Blue to Yellow), a linear path through the RGB cube often intersects the achromatic (gray) axis.
*   **Consequence:** Vibrant input colors result in a desaturated, "muddy" transition in the middle of the gradient.

### B. Mach Banding
Standard multi-stop gradients function as polylines ($C^0$ continuity). At each color stop, the rate of change (velocity) of the color shifts instantaneously.
*   **Consequence:** The human visual system is highly sensitive to discontinuities in the first derivative. These abrupt changes in slope are perceived as physical lines or bands, known as **Mach Bands**, even when the color values are continuous.

### C. Perceptual Non-Uniformity
The sRGB color space is not perceptually uniform. A Euclidean distance of $X$ in the Green channel represents a significantly different visual magnitude than the same distance in the Blue channel.
*   **Consequence:** Gradients appear to accelerate and decelerate unevenly, creating an unbalanced visual rhythm.

---

## 2. Methodology

This tool addresses these limitations by combining a perceptually uniform color space with higher-order interpolation logic.

### The Space: OKLCH
Calculations are performed in **Oklch** (Oklab Cylindrical). By separating Lightness and Chroma from Hue, we can interpolate rotationally. This preserves the perceived intensity and saturation of the colors throughout the transition.

### The Path: Cubic Hermite Splines
Instead of connecting stops with straight lines, we fit a **Catmull-Rom Spline** through the control points.

1.  **$C^1$ Continuity:** The spline ensures that the incoming velocity vector at a stop matches the outgoing velocity vector.
    *   *Result:* Smooth tangents eliminate sharp derivative changes, preventing Mach banding artifacts.
2.  **Curvature:** The spline path naturally curves around the achromatic center rather than cutting through it.
    *   *Result:* Saturation is preserved. For example, a Blue-to-Yellow gradient follows an arc through valid chromatic space, producing a vibrant intermediate tone rather than gray.

> **Note on Optimization:** While this approach is not a global energy minimization curve (Euler-Bernoulli elastica), the Cubic Hermite Spline is a standard local approximation in computer graphics. Since visual smoothness is primarily determined by local continuity between neighbors, this method yields visually superior results without the computational overhead of solving complex differential equations.

---

## 3. The Visual Editor

Editing 3D coordinates on a 2D screen is inherently difficult. The web-based editor solves this via **Smart Projection**.

*   **Osculating Plane Projection:** When a stop is selected, the viewport automatically aligns to the local curvature plane defined by the active point and its neighbors. This ensures edits are made in the most relevant geometric context.
*   **Adaptive Scaling:** The editor uses an $L_2$ norm scaling factor to ensure the visualization remains stable and usable regardless of the curve's rotation in 3D space.

### Tangent Control Modes
*   **Free Mode:** Allows manipulation of both the direction and magnitude (tension) of the tangent vector.
*   **Strength Mode:** Constrains manipulation to the magnitude only. This allows for adjusting the "tension" of the curve—tightening or loosening the transition—without altering the hue trajectory.

---

## 4. Usage in Production

Native browser support for spline-based interpolation is currently unavailable. To utilize these gradients in production, the curve is sampled at discrete intervals (e.g., 40 steps) to generate a standard CSS `linear-gradient` string.

Because OKLCH interpolation in CSS is linear between stops, utilizing a sufficient number of steps renders the linear segments visually indistinguishable from the ideal mathematical curve.

```css
/* Output is a high-resolution linear approximation */
background: linear-gradient(to right, 
  oklch(0.6 0.15 240) 0%, 
  oklch(0.62 0.16 235) 2.5%,
  /* ... intermediate samples ... */
  oklch(0.8 0.1 100) 100%
);
```

## 5. Standalone Library

The repository includes a standalone JavaScript module (`spline.js`) for programmatic generation.

**Note:** The standalone library supports configuring tangent **strengths** (tension) but abstracts away arbitrary 3D direction vectors. In practice, desired path alterations are best achieved by inserting intermediate color stops rather than manually manipulating 3D tangent vectors.
