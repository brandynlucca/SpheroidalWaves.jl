---
title: 'SpheroidalWaves.jl: Arbitrary-precision prolate and oblate spheroidal wave functions in Julia.'
tags:
  - Julia
  - Fortran
  - spheroidal wave functions
  - numerical analysis
  - scientific computing
  - arbitrary precision
authors:
  - name: Brandyn Lucca
    orcid: 0000-0003-3145-2969
    affiliation: '1'
affiliations:
  - name: Fisheries Collaborative Program, University of California Santa Cruz, 1156 High Street, Santa Cruz, CA 95064
    index: 1
date: 17 November 2026
bibliography: paper.bib
---

# Summary

[`SpheroidalWaves.jl`](https://github.com/brandynlucca/SpheroidalWaves.jl) is a Julia package for computing prolate and oblate spheroidal wave functions [@lucca_spheroidalwaves]. Spheroidal wave functions arise when solving the Helmholtz equation by separation of variables in prolate or oblate spheroidal coordinates [@dlmf{}, §§[30.13(iv)](https://dlmf.nist.gov/30.13#iv) and [30.14(iv)](https://dlmf.nist.gov/30.14#iv)]. The package evaluates angular and radial functions for real and complex arguments and provides eigenvalues, expansion coefficients, selectable normalizations, scaled outputs, and coordinate derivatives through fourth order. Double, quadruple, and arbitrary-precision arithmetic support numerical reference calculations. Concentration eigenvalues also quantify the fraction of energy retained within a finite interval by band-limited modes.

# Statement of need

Numerical analysts and researchers solving differential equations or studying spectral concentration need reliable spheroidal wave functions across orders, parameters, and arithmetic precisions. Cancellation and truncation errors require separate refinement of arithmetic precision and expansion length. Established Fortran solvers [@vanburen_mswf] and arbitrary-precision C++ software [@adelman2014] address function evaluation, but projection, parameter differentiation, and integral-operator calculations also require compatible coefficients and normalization conventions. `SpheroidalWaves.jl` supplies function and derivative calculations across real and complex arguments through a reusable Julia interface with shared normalization and precision controls. Angular regularity fixes the shared separation eigenvalue $\lambda_{mn}(c)$ for order $m$, degree $n$, and spheroidal parameter $c$ [@dlmf]. Angular and radial functions $S_{mn}(c,\eta)$ and $R_{mn}(c,x)$ use coordinates $\eta$ and $x$, respectively. The [mathematical documentation](https://brandynlucca.github.io/SpheroidalWaves.jl/dev/math-and-usage/) gives governing equations, normalization conventions, and coordinate domains for the dimensionless quantities.

# State of the field

The [`MathieuandSpheroidalWaveFunctions`](https://github.com/MathieuandSpheroidalWaveFunctions) collection (hereby abbreviated as `MSWF`) supplies the Fortran numerical foundation [@vanburen_mswf]. Prolate methods evaluate first-kind [@vanburen_boisvert2002] and second-kind radial functions and derivatives [@vanburen_boisvert2004]. Oblate methods cover real [@vanburen2017] and complex arguments [@vanburen2020]. The inherited solvers provide degree and angular batching, accuracy estimates, and compiler-selected double, quadruple, and hybrid arithmetic.

Precision and argument domains constrain available reference calculations (Table \ref{tab:functions}). [`SpheroidalWaveFunctions.jl`](https://github.com/uvegege/SpheroidalWaveFunctions.jl) (hereby abbreviated as `SWFuncs.jl`) uses `Float64` [@uvegege_swf], and [`SciPy`](https://docs.scipy.org/doc/scipy/reference/special.html#spheroidal-wave-functions) supplies real-parameter functions [@2020SciPy-NMeth]. [`echoSMs`](https://github.com/ices-tools-dev/echoSMs) [@echosms] and [`acousticTS`](https://github.com/brandynlucca/acousticTS) [@lucca_acousticts] expose prolate calculations in Python and R. The [`scattering`](https://github.com/radelman/scattering) library provides C++/MPFR arbitrary precision and MATLAB interfaces [@adelman2014]. [\texttt{Wolfram Language}](https://reference.wolfram.com/language/guide/SpheroidalAndRelatedFunctions.html) supports arbitrary precision and complex arguments [@wolfram_mathematica], with related methods described by @falloon2003. [`SpheroidalWaves.jl`](https://github.com/brandynlucca/SpheroidalWaves.jl) combines complex-argument and arbitrary-precision evaluation with coefficient access and derived calculations for Julia applications.

`SpheroidalWaves.jl` combines function evaluation, coefficients, derivatives, and integral eigenvalues under shared normalization and precision controls. Reusing `MSWF` retains established real- and complex-parameter solvers. Julia routines add coefficient access, parameter derivatives, and integral eigenvalues to the Fortran evaluations. Extending `SWFuncs.jl` would require generalizing `Float64` calculations, while `scattering` would require complex-parameter algorithms to support the same domains.

```{=latex}
\begin{table}[ht]
    \centering
    \caption{
        Function coverage and supported arithmetic. Superscripts on $S$ and $R$ denote kinds. D, Q, and A denote double, quadruple, and arbitrary precision, respectively. All packages support prolate functions and real $c$. Throughout, $\checkmark$ denotes a supplied function or qualified general operator, and $\times$ denotes no dedicated interface.
    }
    \label{tab:functions}
    \small
    \setlength{\tabcolsep}{3pt}
    \renewcommand{\arraystretch}{1.15}
    \resizebox{\linewidth}{!}{%
        \begin{tabular}{|*{1}{r|}*{1}{l|}*{7}{c|}}
            \hline
                \textbf{Software} & 
                \textbf{Language(s)} & 
                $\boldsymbol{S}_{\mathrm{mn}}^{(1)}$ & 
                $\boldsymbol{S}_{\mathrm{mn}}^{(2)}$ & 
                $\boldsymbol{R}_{\mathrm{mn}}^{(1,2)}$ & 
                $\boldsymbol{R}_{\mathrm{mn}}^{(3,4)}$ & 
                \textbf{Oblate} & 
                \textbf{Complex $c$} & 
                \textbf{Precision} \\
            \hline
                \texttt{echoSMs} & $\text{Python}^a$ & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & $\times$ & $\times$ & D \\
            \hline
                \texttt{SWFunctions.jl} & Julia & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & D \\
            \hline
                \texttt{SciPy} & Python, C++ & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & D \\         
            \hline
                \texttt{acousticTS} & R, C++, Fortran & $\checkmark$ & $\times$ & $\checkmark$ & $\checkmark$ & $\times$ & $\times$ & D, $\text{Q}^b$ \\ 
            \hline
                \texttt{MSWF} & Fortran & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & $\checkmark$ & $\checkmark$ & D, $\text{Q}^c$ \\ 
            \hline
                \texttt{scattering} & C++, MATLAB & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & $\checkmark$ & $\times$ & $\text{A}^d$ \\ 
            \hline
                \texttt{Wolfram} & Wolfram Language & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\times$ & $\checkmark$ & $\checkmark$ & D, A \\ 
            \hline
                \texttt{SpheroidalWaves.jl} & Julia, Fortran, C & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & D, Q, A \\ 
            \hline
        \end{tabular}
    }
    \vspace{4pt}
    \begin{minipage}{\linewidth}
        \fontsize{8}{9.5}\selectfont
        \textsuperscript{a}Python binding to \texttt{MSWF}, maintained by the \texttt{echoSMs} authors. \\
        \textsuperscript{b}Optional \texttt{binary128} computation with \texttt{binary64} R inputs and outputs. \\
        \textsuperscript{c}\texttt{MSWF} defaults to \texttt{selected\_real\_kind(8)}, conventionally double precision. Configurable compiler-supported kinds and hybrid calculations do not provide arbitrary-bit arithmetic. \\
        \textsuperscript{d}Arbitrary-precision kernels and decimal-text output, with \texttt{binary64} MATLAB loaders.
    \end{minipage}
\end{table}
```

Sensitivity analysis requires derivatives with respect to the spheroidal parameter alongside coordinate derivatives (Table \ref{tab:derivatives}). All compared packages supply first coordinate derivatives. `SciPy` and `Wolfram` provide general differentiation operators, while `SpheroidalWaves.jl` supplies dedicated parameter sensitivities under shared normalization and precision controls. Consistent conventions allow function values and derivatives to enter the same calculation without separate rescaling.

```{=latex}
\begin{table}[H]
    \centering
    \caption{
        Derivative coverage for angular or radial solution $f$, with primes denoting coordinate derivatives and $\mathbf J$ denoting a Jacobian. Spheroidal wave Jacobians include $\partial_c f$ and $\partial_c f'$. Checks cover real $c$, with separate real- and imaginary-direction partials available in \texttt{Wolfram} and \texttt{SpheroidalWaves.jl} on smooth complex branches.
    }
    \label{tab:derivatives}
    \small
    \resizebox{\linewidth}{!}{%
        \begin{tabular}{|r|*{8}{l|}}
            \hline
                \textbf{Metric} & 
                \rotatebox{90}{\hspace{4pt}\texttt{MSWF}\hspace{4pt}} & 
                \rotatebox{90}{\hspace{4pt}\texttt{SWFunctions.jl}\hspace{4pt}} & 
                \rotatebox{90}{\hspace{4pt}\texttt{SciPy}\hspace{4pt}} & 
                \rotatebox{90}{\hspace{4pt}\texttt{echoSMs}\hspace{4pt}} & 
                \rotatebox{90}{\hspace{4pt}\texttt{acousticTS}\hspace{4pt}} & 
                \rotatebox{90}{\hspace{4pt}\texttt{spheroidal}\hspace{4pt}} & 
                \rotatebox{90}{\hspace{4pt}\texttt{Wolfram}\hspace{4pt}} & 
                \rotatebox{90}{\hspace{4pt}\texttt{SpheroidalWaves.jl}\hspace{4pt}} \\
            \hline
                Coordinate $f^\prime$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ & $\checkmark$ \\
            \hline
                Coordinate $f^{\prime \prime}$ & $\times$ & $\times$ & $\checkmark^a$ & $\times$ & $\times$ & $\times$ & $\checkmark^b$ & $\checkmark$ \\
            \hline
                $\mathbf{J}_{(\eta,\,c)}(S_{mn})$ & $\times$ & $\times$ & $\checkmark^a$ & $\times$ & $\times$ & $\times$ & $\checkmark^a$ & $\checkmark$ \\
            \hline
                $\mathbf{J}_{(x,\,c)}(R_{mn})$ & $\times$ & $\times$ & $\checkmark^a$ & $\times$ & $\times$ & $\times$ & $\checkmark^a$ & $\checkmark$ \\ 
            \hline
                $\mathbf{J}_{c}(\lambda_{mn})$ & $\times$ & $\times$ & $\checkmark^a$ & $\times$ & $\times$ & $\times$ & $\checkmark^{a,\,b}$ & $\checkmark$ \\
            \hline
        \end{tabular}
    }
\end{table}
\nopagebreak[4]
\vspace{-8pt}
\begingroup
\fontsize{8}{9.5}\selectfont
```

\textsuperscript{a}General numeric differentiation in [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.differentiate.derivative.html) and [Wolfram](https://reference.wolfram.com/language/NumericalCalculus/ref/ND.html).\
\textsuperscript{b}Symbolic differentiation includes [separation-eigenvalue derivatives](https://reference.wolfram.com/language/ref/SpheroidalEigenvalue.html).

```{=latex}
\par
\endgroup
```

Arithmetic precision does not establish accuracy after cancellation and expansion truncation (Table \ref{tab:diagnostics}). `MSWF` returns accuracy estimates, and `scattering` checks Wronskians and alternative representations. `SpheroidalWaves.jl` reports separate value and derivative estimates from refinement, including unavailable or failed-refinement status. The diagnostics support convergence assessment without certifying error bounds.

```{=latex}
\begin{table}[H]
    \centering
    \caption{
        Consistency diagnostics, including the radial Wronskian $\mathcal{W}_R=R_{mn}^{(1)}(R_{mn}^{(2)})'-(R_{mn}^{(1)})'R_{mn}^{(2)}$. \texttt{SpheroidalWaves.jl} supplies normalized Wronskian residuals.
    }
    \label{tab:diagnostics}
    \small
    \resizebox{\linewidth}{!}{%
        \begin{tabular}{|r|*{8}{l|}}
            \hline
                \textbf{Diagnostic} &
                \rotatebox{90}{\hspace{4pt}\texttt{MSWF}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{SWFunctions.jl}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{SciPy}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{echoSMs}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{acousticTS}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{spheroidal}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{Wolfram}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{SpheroidalWaves.jl}\hspace{4pt}} \\
            \hline
                $\mathcal{W}\left(R_{mn} \right)$ & $\times^{a}$ & $\times$ & $\times$ & $\times$ & $\times$ & $\checkmark^{b}$ & $\checkmark^{c}$ & $\checkmark$ \\
            \hline
                Accuracy diagnostics & $\checkmark$ & $\times$ & $\times$ & $\times^{d}$ & $\times^{d}$ & $\checkmark^{b}$ & $\checkmark^{e}$ & $\checkmark^{f}$ \\
            \hline
        \end{tabular}
    }
\end{table}
\nopagebreak[4]
\vspace{-8pt}
\begingroup
\fontsize{8}{9.5}\selectfont
```

\textsuperscript{a}Internal calculation.\
\textsuperscript{b}Wronskian residuals and differences between alternative representations.\
\textsuperscript{c}General [Wronskian operator](https://reference.wolfram.com/language/ref/Wronskian.html).\
\textsuperscript{d}Backend accuracy arrays omitted from returned results.\
\textsuperscript{e}General [precision tracking](https://reference.wolfram.com/language/ref/Accuracy.html).\
\textsuperscript{f}Refinement-based estimates.

```{=latex}
\par
\endgroup
```

Separation eigenvalues determine the differential-equation solutions, whereas concentration eigenvalues quantify finite-interval energy retention [@slepian_pollak1961]. `SpheroidalWaves.jl` computes concentration and finite Fourier eigenvalues for real prolate $m=0$, $c\geq0$ (Table \ref{tab:spectral}). `SciPy` and `Wolfram` supply general root finders, while `SpheroidalWaves.jl` provides dedicated eigenvalue inversion. Logarithmic and complementary concentration targets preserve small values near zero and unit concentration.

```{=latex}
\begin{table}[H]
    \centering
    \caption{
        Spectral interfaces, with concentration eigenvalue $\Lambda_n$ and finite Fourier eigenvalue $\mu_n$. \texttt{SpheroidalWaves.jl} inversion requires a real bracket and fixed $(m,n)$.
    }
    \label{tab:spectral}
    \small
    \resizebox{\linewidth}{!}{%
        \begin{tabular}{|r|*{8}{l|}}
            \hline
                \textbf{Spectral calculation} &
                \rotatebox{90}{\hspace{4pt}\texttt{MSWF}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{SWFunctions.jl}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{SciPy}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{echoSMs}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{acousticTS}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{spheroidal}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{Wolfram}\hspace{4pt}} &
                \rotatebox{90}{\hspace{4pt}\texttt{SpheroidalWaves.jl}\hspace{4pt}} \\
            \hline
                Separation $\lambda_{mn}$ & $\times^{a}$ & $\checkmark$ & $\checkmark$ & $\times$ & $\times$ & $\checkmark$ & $\checkmark$ & $\checkmark$ \\
            \hline
                Concentration $\Lambda_n$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\checkmark$ \\
            \hline
                Fourier $\mu_n$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\checkmark$ \\
            \hline
                Inversion $\lambda_{mn}\mapsto c$ & $\times$ & $\times$ & $\checkmark^{b}$ & $\times$ & $\times$ & $\times$ & $\checkmark^{b}$ & $\checkmark$ \\
            \hline
                Eigenmode continuation & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\times$ & $\checkmark$ \\
            \hline
        \end{tabular}
    }
\end{table}
\nopagebreak[4]
\vspace{-8pt}
\begingroup
\fontsize{8}{9.5}\selectfont
```

\textsuperscript{a}Internal calculation and diagnostic-file output.\
\textsuperscript{b}General root finding in [SciPy](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.root_scalar.html) and [Wolfram](https://reference.wolfram.com/language/ref/FindRoot.html).

```{=latex}
\par
\endgroup
```

# Software design

The functions `smn` and `rmn` evaluate angular and radial functions, while `spheroidal` combines wave and coefficient results. Fortran supplies double and quadruple arithmetic, with `MPFR` [@fousse2007] and `MPC` [@enge2022] supporting arbitrary precision. Explicit precision requests refine arithmetic and expansion length separately and preserve requested `BigFloat` precision independently of ambient settings. Scaled mantissa/exponent outputs extend range without increasing accuracy. Unresolved refinement raises an error.

Expansion coefficients give the weights used to reconstruct first-kind angular functions from Ferrers functions. Default normalization and basis choices give [@dlmf{}, §§[30.8(i)](https://dlmf.nist.gov/30.8#i), [30.11(i)](https://dlmf.nist.gov/30.11#i), and [30.11(v)](https://dlmf.nist.gov/30.11#v)]:
$$
  \begin{aligned}
    S_{mn}^{(1)}(c,\eta)
      &=\sum_{\ell\ge m}'d_\ell^{mn}(c)\,\mathsf P_\ell^m(\eta),\\
    A_{mn}(c)&=\sum_{\ell\ge -m}'d_\ell^{mn}(c),\\
    R_{mn}^{(1)}(c,x)
      &=K_{mn}(c)\sum_{\ell\ge m}'d_\ell^{mn}(c)\,P_\ell^m(z),\qquad c\ne0,
  \end{aligned}
  \label{eq:coefficients}
$$
where $d_\ell^{mn}$ (`dmn`) weights each Ferrers function, $A_{mn}$ (`amn`) is a coefficient sum used in angular–radial connection formulas, and $K_{mn}$ (`kmn`) scales the exterior angular expansion into the first-kind radial function. All three quantities are dimensionless. Primed sums select $\ell$ with the parity of $n$. The functions $\mathsf P_\ell^m$ and $P_\ell^m$ are Ferrers and exterior Legendre functions, with $z=x$ (prolate) or $z=-ix$ (oblate). The $A_{mn}$ sum includes coefficients extended below degree $m$ by recurrence. `dmn` also provides parameter derivatives and an orthonormal basis independently of wave normalization. Angular first-kind normalization is Meixner–Schäfke or unit norm for real parameters, analytically continued for complex parameters. Radial kinds one and two offer standard normalization or static normalization with finite $c=0$ limits. `eigenvalue` returns separation, concentration, or Fourier eigenvalues, reusing angular coefficients for the integral operators [@osipov2013]. `find_c_for_eigenvalue` inverts separation or concentration eigenvalues within a real bracket.

Derivative calculations reuse analytic expressions to limit subtraction errors. First parameter derivatives follow coefficient expansions and governing equations. Second wave parameter derivatives difference analytic first derivatives, requiring additional evaluations and step-size refinement. Coordinate derivatives through fourth order follow differentiated equations, with one-sided endpoint limits. Complex-parameter continuation matches angular profiles and phase across derivative samples and reports unresolved continuation. `angular_zeros` locates first-kind angular zeros and stationary points for real $c$. `radial_zeros` locates zeros and stationary points of radial kinds one and two for real nonzero $c$ on finite nonsingular intervals.

Input checks enforce integer order and degree, valid coordinates, and supported option combinations. Requested precision selects `Float64` or `BigFloat` for real and imaginary parts of returned values. Thread-safe evaluation keeps working data separate between simultaneous calculations, at the cost of additional memory. Tests compare simultaneous and sequential calls and verify preservation of unrelated `BigFloat` precision and rounding. [`SpheroidalWaves_jll`](https://github.com/JuliaBinaryWrappers/SpheroidalWaves_jll.jl) and [`Yggdrasil`](https://github.com/JuliaPackaging/Yggdrasil/blob/master/S/SpheroidalWaves/build_tarballs.jl) distribute compiled libraries, reducing compiler requirements.

# Research impact statement

The independently developed `SphericalScattering.jl` package [@hofmann2023] depends on `SpheroidalWaves.jl` for acoustic scattering from prolate and oblate spheroids. The documented [modal series solution for spheroids](https://hobezwe.github.io/SphericalScattering.jl/stable/numerics/spheroidSeries/) uses angular and radial functions for incident-field expansion and sound-hard or sound-soft boundary matching. `AcousticScattering.jl` also integrates `SpheroidalWaves.jl` into modal series solution for acoustic scattering problems [@lucca_acousticscattering], including [fluid-filled spheroids](https://brandynlucca.github.io/AcousticScattering.jl/stable/models/modal/) with coupling between angular degrees. Core elements of the Fortran kernels were also adapted into the `acousticTS` package's methods directly.

Test cases used `SpheroidalWaves.jl` to validate the outputs of the [spheroidal wave functions implemented in `acousticTS`](https://brandynlucca.github.io/acousticTS/reference/index.html#spheroidal-wave-functions) [@lucca_acousticts]. `SpheroidalWaves.jl` was also used for diagnosing angular and radial function calculations from the `Prol_Spheroid` Python library [@khodabandeloo2025] as well as those from [`echoSMs`](https://ices-tools-dev.github.io/echoSMs/api_utilities/#echosms.utils.pro_ang1) and `scattering`. The comparisons assess downstream predictions, while shared `MSWF` methods limit independence at the function-kernel level. Tests compare angular and radial values with stored `Wolfram|Alpha` references and radial values and derivatives with published `MSWF` sample tables [@vanburen_mswf]. Independently generated high-precision expansions, exact limits, and branch checks address errors potentially hidden by Wronskian agreement. The combined evidence supports research use without establishing accuracy throughout the parameter domain.

# AI usage disclosure

AI tools were used in a limited capacity during later package development. OpenAI Codex (GPT-6) was assisted with debugging and refining the `Yggdrasil` workflow ([`Yggdrasil.yml`](https://github.com/brandynlucca/SpheroidalWaves.jl/blob/v0.7.0/.github/workflows/Yggdrasil.yml)) for generating the precompiled binaries. It also assisted with addressing bugs in the continuous integration workflow ([`CI.yml`](https://github.com/brandynlucca/SpheroidalWaves.jl/blob/v0.7.0/.github/workflows/CI.yml)) that arose due to differences in Fortran compilers, precision arithmetic, and parallelization across Windows, Linux, and MacOS containers. All AI-generated test results, artifacts, and reports were systematically human-reviewed and -validated.

# Acknowledgements

The original Fortran-based methods were originally developed by Arnie Lee Van Buren and Jeffrey E. Boisvert, which make up the `MSWF` collection. Their open-source solvers make the recurrence and expansion methods available for reuse in research software, including serving as  the basis for the Fortran kernels in `SpheroidalWaves.jl`. This research received no funding, and the author declares no competing interests.

# References
