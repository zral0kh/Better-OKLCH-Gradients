# Better OKLCH Spline Gradient

We provide both an js implementation and a web-based editor for creating smooth color gradients.
The webapp allows for generating perceptually smooth, high-order color gradients using **Cubic Hermite Splines (Catmull-Rom)** in the **OKLCH** color space.

View the interactive Editor [here](https://zral0kh.github.io/hermite-oklch-gpage/)

## 1. The Problem with Standard Gradients

Most web gradients (CSS `linear-gradient`) rely on **Linear Interpolation (Lerp)** in the **sRGB** color space. While computationally cheap, this approach has three fundamental flaws:

### A. The "Gray Dead Zone"
When interpolating between two complementary colors (e.g., Blue to Yellow), linear interpolation draws a straight line through the 3D color cube. This line often passes directly through the desaturated gray center of the color space.
*   **Result:** Vibrant colors turn into "mud" in the middle of the transition.

### B. Mach Bands (Visual Hard Edges)
Standard multi-stop gradients are "Polylines." At every color stop, the rate of color change (velocity) changes instantly.
*   **Math:** They have $C^0$ continuity (connected) but not $C^1$ continuity (smooth tangent).
*   **Result:** The human eye is incredibly sensitive to changes in the rate of change (derivative). We perceive these sharp velocity changes as physical lines or bands, known as **Mach Bands**, even if the color values are technically continuous.

### C. Perceptual Non-Uniformity
sRGB is not perceptually uniform. A mathematical change of `+10` in Green looks much more drastic to the human eye than `+10` in Blue.
*   **Result:** Gradients appear to speed up and slow down unevenly, making the transition look unbalanced.

---

## 2. The Solution: Oklch + Splines

This tool solves these problems by combining two advanced techniques:

### A. The Space: OKLCH
We perform all calculations in **Oklch** (Oklab Cylindrical).
*   **Perceptual Uniformity:** A distance of $X$ in Oklch space represents the same amount of visual change regardless of the color.
*   **Hue Preservation:** By separating Chroma (Saturation) and Hue, we can interpolate rotationally.
*   **Result:** Smooth, predictable transitions that match how the human eye perceives color.

This is how many people choose the gradients currently, but this still suffers from the issues of linear interpolation mentioned above.


### B. The Path: Cubic Hermite Splines
Instead of drawing straight lines between stops, we fit a **Catmull-Rom Spline** through the points in 3D space.

1.  **$C^1$ Continuity:** The spline ensures that the incoming velocity at a stop perfectly matches the outgoing velocity. There are no sharp corners.
    *   *Benefit:* **Mach Bands are eliminated.** The gradient flows like a liquid rather than a series of hard turns.
2.  **Curved Paths:** A spline can curve *around* the gray center of the color space rather than cutting through it.
    *   *Benefit:* **Preserved Saturation.** A gradient from Blue to Yellow can curve slightly outward to maintain vibrancy, resulting in a rich intermediate green or white/pink (depending on the path) instead of gray mud.

In general this is not globally optimal for larger than 4-stop gradients, because it does not minimize total curvature, only locally.
Nontheless, because humans are more sensitive to local changes, this still produces a visually superior result in practice.

---

## 3. Key Features & Math

### 3D Tangent Control
In this tool, a color stop is not just a point; it is a point with a **velocity vector (tangent)**.
*   **Direction:** Controls where the gradient goes next (e.g., "head towards Cyan before turning to Blue").
*   **Magnitude (Tension):** Controls how "aggressively" the color moves. A long tangent makes the color "linger" longer before accelerating to the next stop.

### Smart Projection (Osculating Plane)
Editing 3D coordinates on a 2D screen is notoriously difficult.
*   **The Solution:** When you select a stop, the editor projects the view onto the **local curvature plane** (defined by the active stop and its neighbors).
*   **Why it matters:** This allows you to edit the curve naturally without accidental distortions. The $L_2$ norm scaling ensures the zoom level remains stable regardless of the camera's rotation.

### Editor Modes
*   **Free Mode:** Allows changing both the direction and the magnitude of the color velocity.
*   **Strength (Locked) Mode:** Constrains changes to the magnitude only. This effectively acts as a **Tension** control, allowing you to tighten or loosen the curve without altering the hue path.

---

## 5. Usage in Production

Browsers do not yet support `spline-gradient()` natively. To use these gradients on the web, we sample the polynomial curve at discrete intervals (e.g., 40 steps) and output a standard linear-gradient string.

```css
/* Resulting CSS approximates the curve using many linear stops */
background: linear-gradient(to right, 
  oklch(0.6 0.15 240) 0%, 
  oklch(0.62 0.16 235) 2.5%,
  /* ... intermediate steps ... */
  oklch(0.8 0.1 100) 100%
);
```

Because Oklch interpolation in CSS is linear between stops, using enough steps makes the linear segments indistinguishable from the true curve to the human eye.

## Try it Out
You may either use the [online editor here](https://zral0kh.github.io/hermite-oklch-gpage/) or download the js file and use it in your own projects!
The js file currently only allows setting strengths and not directions of the tangents. You can just use different color points instead to achieve the same result in theory, may need some tweaking. In most cases it should not be necessary to change tangent directions anyway.
