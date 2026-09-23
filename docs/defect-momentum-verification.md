# Single-defect momentum spectra: completed verification

The 1D, zero-temperature, one-electron calculation uses t = omega0 = U = 1,
lambda = 0.25 and 0.5, and a normalized coherent probe on 33 sites centered on
the defect. All 201 momenta in [0, pi] are directly contracted from nonlocal
propagation or sampled endpoint displacements. The 1601 displayed energies are
not an assertion of inverse-Laplace resolution.

- VED was completed first, with clean/defect maps, an additional 65-site probe,
  independent cloud/radius/recursion checks and direct cosine/sine injections.
- Fresh nonlocal bare DiagMC used 16 independent chains and 2.048 billion production
  steps on highmem. No order cap was reached. Three analytic limits passed before
  production; their maximum standardized discrepancy was 2.211.
- All 402 reference-free continuation selections and 12864 conditional bootstrap
  repeats completed. The largest selected-fit moment residual was 9.295e-8.
- Uncontinued Green functions differ from VED by at most 2.471 conservative
  standard errors for tau <= 6. At eta = 0.25, full-spectrum relative L1
  differences are 3.97–12.13% and 6.23–15.67% for the two couplings. At eta = 1,
  the respective maxima are 1.685% and 1.749%. Narrow features remain unvalidated.
- Fresh extraction verified 2215 file hashes across 24 archives and 27 direct
  assets. Compiled update tests and all 26 Python tests passed. Six representative
  spectra and their bootstrap quantiles were regenerated, and an analytic
  fixed-seed sampler replay produced bitwise-identical raw blocks.
- Website build/check passed for eight bilingual pages, 700 LaTeX expressions
  with MathML, four existing TikZ diagrams and local scientific resources.
  Browser checks passed at widths 1280, 390 and 320 with JavaScript disabled,
  external requests blocked, correct language switching and no document overflow.
  The new section's four figures also passed focused checks in both languages.

The pages embed PNG previews to keep the Pages artifact compact; original vector
SVG/PDF and NPZ data are in the [release](https://github.com/rjguo1208/Holstein-model/releases/tag/defect-momentum-20260923).
The [fresh-extraction record](https://github.com/rjguo1208/Holstein-model/releases/download/defect-momentum-20260923/verification.json)
and [file manifest](https://github.com/rjguo1208/Holstein-model/releases/download/defect-momentum-20260923/manifest.json)
provide reproducibility details. Numerical methods and resolution qualifications
are documented in [DEFECT_MOMENTUM.md](../research/holstein-defect/DEFECT_MOMENTUM.md).

中文：本轮先完成201点 VED 动量谱，再独立完成20.48亿步非局域 bare DiagMC
重新采样及12864次条件 bootstrap。2215个文件的哈希、全新解压后的编译与26项
Python 测试、代表性谱及固定种子轨迹复现均通过。中英文页面在桌面和手机宽度下
通过公式、图片、语言切换和页面溢出检查。图像的密集网格不等于细峰已经分辨；
窗口、声子截断、统计噪声及延拓误差的限制均保留在网页和方法说明中。
