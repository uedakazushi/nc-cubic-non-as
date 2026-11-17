# 量子Fermat三次曲面はAS surfaceではない — nc-cubic-non-as

**最初に [研究方針](00_research_plan.md) を読む。** 公表用の独立した日本語LuaLaTeX短報は [quantum_fermat_not_as_surface.tex](papers/quantum_fermat_not_as_surface.tex)、現在の研究ノートは [research/research_notes.md](research/research_notes.md)。

量子Fermat三次曲面の導来圏を非自明な離散torsionを持つ有限群商stackの捩れ付き層の導来圏と比較し、固定楕円曲線のtwisted orbifold HKRから `dim HH_{±1}=2`, `dim HH_0=5` を示す。OUのAS surfaceには完全強例外列があるため、この曲面はその定義のAS surfaceではない。**同論文の証明と検算の範囲を区別し、全H3の分類や新規性を完了済みとはしない。**

## 再現

```sh
python3 -m pip install -r requirements.txt
python3 computations/verify_all.py
latexmk -lualatex papers/quantum_fermat_not_as_surface.tex
```

GitHub mainが唯一の正本。PDF・ZIP・中間生成物は追跡しない。より広い分類の正本は [nc-cubic-moduli](https://github.com/uedakazushi/nc-cubic-moduli)、比較は [nc-cubic-ambient-moduli](https://github.com/uedakazushi/nc-cubic-ambient-moduli)、楕円E6量子化は [nc-cubic-e6](https://github.com/uedakazushi/nc-cubic-e6)。研究資料は [固定source commit](https://github.com/uedakazushi/nc-cubic-moduli/tree/c5fd8a8301f6b1c283a6f0510cde80ae6a006e80) から移行した。
