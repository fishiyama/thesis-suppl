D論補足説明用のデモプログラム Juliaで作成
https://tsukuba.repo.nii.ac.jp/records/2000757

tutorial.jl	 プログラム本体
tutorial.jl.pdf	 プログラムのプリントアウト
tutorial.out.pdf プログラムの実行結果に解説を記入したもの
プログラムは上記解説記述後の修正が含まれるので、同一にはなっていない。

高階差分方程式を解く方式なので、高階差分が桁落ちしないよう、
隣接サンプル間の差分が大きくなるように注意する必要がある。
特に実測データの解析時には有効桁数の確保が困難なので注意。

ishiyama.fumihiko.xd@alumni.tsukuba.ac.jp
2022/09/27
