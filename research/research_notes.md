# 量子Fermatの非AS障害：現行研究ノート

2026-11-17 / ou-obstruction-v1。研究方針は `00_research_plan.md`、全文証明は `papers/quantum_fermat_not_as_surface.tex`。

## 幾何学と主張

四変数量子多項式環の中心Fermat三次商のqgrの導来圏を、古典Fermat三次曲面XのH=(mu_3)^2商stackに非自明離散torsionを入れたcoherent層の圏に同定する。twisted orbifold HKRでH-fixed sectorのHodge固有指標を取り、HH_{-1}とHH_1が各2次元、HH_0が5次元となる。OU定義のAS surfaceは完全強例外列を持ち、奇数次HHが零。従って本例はOU AS surfaceではない。

## 前提の検証範囲

代数のgrading・Nakayama・中心性・Hilbert値の有限計算は `computations/verify_quantum_fermat.py`、twisted sectorの列挙は `computations/verify_quantum_fermat_ou_sectors.py` が再現する。 `computations/quantum_fermat_ou_sector_certificate.json` に現在のsector別結果を保持する。別のH3認識や正則性の計算は `verify_boundary_fermat.py` が有限範囲を補助する。これらは導来同値や全局所正則性を単独で保証しない。

原稿には既存のOU定義・例外列と本研究の捩れ付き商の比較を記載。公表前には一次文献からcocycle規約、局所Hochschild sectorへの作用、固定楕円曲線のHodge指標、捩れの非退化性を独立監査する。導来等価の具体的functorと生成元を十分に示す必要がある。新規性の徹底調査は未完成。

## 再現

`python3 computations/verify_all.py` は上記2つの主要計算を実行し、HH次元のJSON結果を検査する。数値例から全成分を外挿しない。LuaLaTeXで `latexmk -lualatex papers/quantum_fermat_not_as_surface.tex` を実行する。元の全H3分類を丸ごとコピーしない。
