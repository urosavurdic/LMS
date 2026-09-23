# LMS adaptive filter

Adaptive noise cancellation with the LMS algorithm, comparing a fixed step
size against two variable-step-size variants where μ(n) is driven by the
instantaneous error. Learning curves are averaged over 100 Monte Carlo runs at
each of three SNRs.

## The idea

Standard LMS uses one fixed step size, which forces a trade-off: large enough
to converge quickly, or small enough to sit close to the optimum once it gets
there — not both. The variants make the step size a function of the error, so
it is large while the filter is still far off and shrinks as the error does.

| Filter | Step size |
|---|---|
| LMS | μ = 1e-5, constant |
| Variant 1 | μ(n) = α·log₁₀(1 + ½·(e(n)/δ)²) |
| Variant 2 | μ(n) = α·log₁₀(1 + ½·\|e(n)·e(n−1)\|/δ²) |

Variant 2 uses the product of successive errors rather than the current one
squared. Because the noise component is uncorrelated between samples, that
product tends toward zero while a genuine tracking error persists — so the
step size responds to the trend instead of to one noisy sample.

## Setup

A sinusoid of period 20 samples in additive white Gaussian noise. The filter
input is that signal delayed by one sample, making this a one-step linear
predictor: a 100-tap filter, 10,000 samples per run, 100 runs averaged, at
SNRs of 10, 20 and 30 dB.

## Running it

Open `lms_adaptive_filter.m` in MATLAB and run it. Needs the Communications
Toolbox for `awgn` and the Signal Processing Toolbox for `freqz`.

It produces, for each SNR: the input, noisy and filtered signals; the spectra
of all three; and the averaged mean-square-error learning curves.

| Path | What it is |
|---|---|
| `lms_adaptive_filter.m` | the simulation |
| `report.docx` | written analysis |
| `slides.pptx` | presentation |

## License

MIT — see [LICENSE](LICENSE).
