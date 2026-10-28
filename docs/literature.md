# Earlier work on hemispheric asymmetry, and where this project sits

This is a working map of the literature the project builds on or argues with. It is not a full review. For that see Norton, Charbonneau & Passos (2014) on observations vs dynamo models, Hathaway (2015) on the solar cycle in general and Petrovay (2020) on cycle prediction.

## 1. The chain everyone uses

The Babcock–Leighton picture links three things:

- **A**: active regions of cycle *n* emerge, tilted (Joy's law), unevenly in the two hemispheres.
- **B**: their trailing flux drifts poleward and builds the polar field. Its north–south imbalance at the minimum is the polar asymmetry.
- **C**: the polar (poloidal) field seeds the toroidal field and so the activity of cycle *n+1*, separately in each hemisphere if the hemispheres are weakly coupled.

Most earlier work tests **B → C**. This project measures **A → B** cycle by cycle and also tests **A → C** directly.

## 2. Observations of the asymmetry

| Work | What it shows | Relevance here |
|---|---|---|
| Norton, Charbonneau & Passos 2014, *Space Sci. Rev.* 186, 251 ([arXiv:1411.7052](https://arxiv.org/abs/1411.7052)) | From RGO data, hemispheric sunspot-area production is at most ~20% asymmetric for cycles 12–23. One hemisphere tends to lead for several cycles. Compares with a 1600-yr dynamo simulation and puts limits on cross-equatorial coupling. | The range of A^n our targets live in. Persistence of the leading hemisphere is the baseline channel 2 has to beat (`rev.t4`). |
| Svalgaard & Kamide 2013, *ApJ* 763, 23 ([arXiv:1207.2077](https://arxiv.org/abs/1207.2077)) | Polar reversals are asymmetric because activity is: the more active hemisphere early in the cycle reverses first. | A qualitative A → B statement about *timing*. We test the *amplitude* at the following minimum, and the reversal epoch t_R is the anchor of our window. |
| Muñoz-Jaramillo et al. 2012, *ApJ* 753, 146 ([arXiv:1303.0345](https://arxiv.org/abs/1303.0345)) | Calibrates a century of Mount Wilson polar faculae counts against WSO and MDI polar flux. | Source of our A_pf for SC15–19. |
| Muñoz-Jaramillo et al. 2013, *ApJL* 767 ([arXiv:1304.3151](https://arxiv.org/abs/1304.3151)) | Polar flux at minimum predicts the next cycle; sunspot area alone does not, unless weighted by tilt. Memory is about one cycle. | One-cycle memory is why the predictor uses cycle n only; WS-E checks whether cycle n−2 adds anything. The tilt result motivates the tilt-weighted X2 in Addendum B. |

## 3. Models: where does the asymmetry come from?

| Work | Mechanism | Relevance here |
|---|---|---|
| Goel & Choudhuri 2009, *Res. Astron. Astrophys.* 9, 115 ([arXiv:0712.3988](https://arxiv.org/abs/0712.3988)) | Randomness in the Babcock–Leighton source leaves one hemisphere with a stronger polar field at the end of a cycle, which then shows up as asymmetry in the next cycle. Modelled with the Surya flux-transport dynamo. | The B → C half of our chain, and the direct ancestor of the `dynamo/` extension. |
| Choudhuri, Chatterjee & Jiang 2007, *PRL* 98, 131103 ([astro-ph/0701527](https://arxiv.org/abs/astro-ph/0701527)) | Feeds observed polar fields into the dynamo at each minimum and predicts cycle amplitude. | The obvious way to drive Surya with *hemispheric* observed polar fields. |
| Chatterjee, Nandy & Choudhuri 2004, *A&A* 427, 1019 ([astro-ph/0405027](https://arxiv.org/abs/astro-ph/0405027)) | Full-sphere Surya; dipolar parity needs efficient cross-equatorial diffusion of the poloidal field. | The model in `dynamo/`. Its hemispheric coupling sets how long an imposed asymmetry can survive. |
| Belucz & Dikpati 2013, *ApJ* | Changing the meridional flow in one hemisphere changes that hemisphere's cycle almost alone. | An alternative to the B–L source as the origin. Our WS-A tests hemispheric transport asymmetry in the SFT. |
| Shetye et al. 2015 ([arXiv:1412.3245](https://arxiv.org/abs/1412.3245)) | Observations and flux-transport-dynamo modelling of N–S asymmetries. | Same family as above. |
| Passos et al. 2014, *A&A* ([arXiv:1309.2186](https://arxiv.org/abs/1309.2186)) | Mean-field α plus B–L source with fluctuations gives grand minima/maxima and hemispheric asymmetry. | The poloidal source and noise of the SURYA runs in `extension/` (v7, v8), as used by Bhowmik (2019). |
| Shukuya & Kusano 2017 ([arXiv:1612.03294](https://arxiv.org/abs/1612.03294)) | A flux-transport dynamo reproduces the hemispheric phase asymmetry without an ad hoc asymmetric source. | Related to the phase lag test `rev.n1` and N2 ("does the first-reversing hemisphere set A_pf"). |
| Schüssler & Cameron 2018, *A&A* 618, A89 ([arXiv:1807.10061](https://arxiv.org/abs/1807.10061)) | Asymmetry is a weakly damped quadrupolar mode kicked by random fluctuations of the poloidal source; explains spectral peaks near 8.5 yr and 30–50 yr. | The main theoretical frame: asymmetry is noise in the source term. Our size decomposition asks *which* regions carry that noise. |
| Nepomnyashchikh et al. 2019 ([arXiv:1903.10707](https://arxiv.org/abs/1903.10707)) | Long-term asymmetry from fluctuations in dynamo parameters, through an equator-symmetric field component. | Same frame, long time scales. |
| Kitchatinov & Khlystova 2021, *ApJ* 919, 36 ([arXiv:2107.00865](https://arxiv.org/abs/2107.00865)) | Dipole–quadrupole superposition; phase locking gives runs of northern or southern asymmetry; polar asymmetry at minima as a precursor. | Predicts persistence of the sign. Our SC14 and SC21 sign changes (`cr.z4`, WS-F) are the cases such models have to explain. |

## 4. Surface flux transport and prediction of the asymmetry

| Work | What it shows | Relevance here |
|---|---|---|
| Jiang, Cameron & Schüssler 2015, *ApJL* 808, L28 ([arXiv:1507.01764](https://arxiv.org/abs/1507.01764)) | The weak cycle 24 traces back to a weak polar field, which in turn traces back to a few regions with wrong orientation or small tilt and to cross-equatorial flux. | Individual regions matter. That motivates our per-region (peak-area) size classes. |
| Hathaway & Upton 2016, *JGR* ([arXiv:1611.05106](https://arxiv.org/abs/1611.05106)) | AFT simulations predict a south-dominant polar field at the end of cycle 24. | A direct comparison point for our registered SC25 test (`rev.n7`). |
| Nagy et al. 2017, *Sol. Phys.* ([arXiv:1712.02185](https://arxiv.org/abs/1712.02185)) | A single "rogue" bipolar region can boost, suppress or even stop later cycles; polar-cap asymmetry in cycle i correlates with activity asymmetry in cycle i+1. | The giant-region class in our SFT budget and experiments (`cr.x4`, `cr.x9`, `rev.t9`) is the observational version of this. |
| Nagy, Lemerle & Charbonneau 2020 ([arXiv:1909.07672](https://arxiv.org/abs/1909.07672)) | Rogue regions and hemispheric asymmetry specifically. | Same. |
| Bhowmik 2019, *A&A* 632, A117, [doi:10.1051/0004-6361/201834425](https://doi.org/10.1051/0004-6361/201834425) ([arXiv:1906.06877](https://arxiv.org/abs/1906.06877)) | Over ~110 yr, the polar-flux asymmetry at minimum is closely tied to the next cycle's hemispheric activity. SFT runs point at tilt randomness as the main source of polar asymmetry. | We reproduce the B → C test on the corrected catalogue (`cr.y1`, Fig. 8) and compare every one of her results with ours in [`extension/reports/REPORT_v5_additions.md`](../extension/reports/REPORT_v5_additions.md); v8 reruns our analyses in her SFT and dynamo. |
| Karak & Miesch 2017, *ApJ* 847, 69 ([arXiv:1706.08933](https://arxiv.org/abs/1706.08933)) | Tilt scatter in a 3-D B–L dynamo gives realistic cycle variability. | Why the SFT ensemble carries 101 tilt-scatter members instead of one deterministic run. |

## 5. What this project adds

Put plainly, and limited to what the analysis actually tests:

1. **A → B as a cycle-by-cycle number.** Earlier work says the more active hemisphere reverses first (Svalgaard & Kamide 2013) or shows A → B in SFT runs (Bhowmik 2019). Here a single observed predictor, the latitude-weighted area imbalance in the ~2 years after the polar reversal, tracks A_pf at the next minimum over cycles 15–24 (r = 0.93; 0.86 under nested leave-one-cycle-out).
2. **A → C without going through the polar field.** The same kind of predictor over a wider window tracks the next cycle's area asymmetry over cycles 12–23 (r = 0.77; 0.60 nested). It is tested against plain persistence, the obvious alternative, and most of it turns out to be persistence (the dominant hemisphere tends to stay dominant), which the extension then studies directly.
3. **Per-region size classes.** The corrected catalogue has unique region ids, so each region is classified once by its peak area. Older row-level classifications count a long-lived large region as many small rows.
4. **Model read out like the data.** The SFT ensemble is sampled at the same ±1.6 yr windows and caps as the observed A_pf, so model and observation are compared like for like.
5. **The two sign changes, SC14 and SC21.** They are treated as tests rather than as outliers.
6. **Persistence and flips of the dominant hemisphere.** The area imbalance has a cycle-to-cycle memory (lag-1 r = +0.57, n = 12), and both sign changes go through a near-balanced cycle. A 1-D SFT has no such memory. A flux-transport dynamo with Passos-type noise (Bhowmik's set-up) carries most of it (+0.45), though with shorter same-sign runs than the Sun (`extension/`).
7. **A candidate origin for the parity factor.** Even cycles switch their dominant hemisphere after maximum and odd cycles do not, so declining-phase emergence gets the A_pf sign right in 9 of 10 cycles without the factor. This is exploratory (`extension/`, Q2).
8. **Statistics built for n ≈ 10.** Every scanned quantity carries a nested value and a family-wise p; a specification curve covers ~2 M definitions; errors in A_pf are propagated; BH/Holm correction is applied across all tests; one prediction (SC25) was registered before the data exist. This kind of protocol is rare in the asymmetry literature, where single correlations over ~10 cycles are often reported without a search penalty.

What it does **not** add yet: an established mechanism for the parity factor (point 7 is exploratory), an SFT that works in even cycles (the project's SFT and Bhowmik's both reproduce A_pf in odd cycles only), or a dynamo that reproduces the Sun's long same-sign runs.

## References not linked above

- Hathaway, D. H. 2015, *Living Rev. Sol. Phys.* 12, 4.
- Petrovay, K. 2020, *Living Rev. Sol. Phys.* 17, 2 ([arXiv:1907.02107](https://arxiv.org/abs/1907.02107)).

Journal details above were checked against abstracts and repository records, not the full papers. Check volume and page numbers before citing them in a manuscript.
