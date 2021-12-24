# 球座標系2次元FDTD法を用いたVLF/LF帯電磁波による下部電離圏電子密度同定

## コード改変の際のローカルルール
* ソースには誰が見ても分かるコメントをつける (日本語で可)
 * なるべくコメントが必要のない変数名・関数名をつけること。a, b, ..とか、t1, t2とかはなるべく避ける
 * コメントは物理量名以外に、単位をつける
 * コメントは /* */ で書くようにする。// のコメントは「このコードを動作させない」意味で用いる
* コミットする際のコメントも誰が見ても分かるようにつけること
* コードを改変したら動作を確認してからコミットすること

## 次のタスク
1. 現在は、f=40kHz、伝搬距離2,200kmの計算をしているので、計算時間が長くかかります。f=22.2kHz、伝搬距離1,000kmにしましょう。fdtd2d.hを書き換えて下さい。また、その際に<b>送信局をえびの(JJI)</b>に直して下さい(緯度経度と磁場)。もう以前に渡したファイルを使って、それで計算をして、地表面電界強度と位相を計算しておいて下さい。
2. クラスの概念がつかめていないようなので、単なる関数でやりましょう。
 * main関数は z', beta, Lp, z_dec, sig_per の「電離圏のパラメタ」をFDTD計算関数に渡し、関数は地表面の電界(22.2kHzのスペクトル)を返してくる形にしたいと思います。従って、main関数は以下のようになります。
```c++:main.cpp
int main(void){
  double z_prime = 85.0e3; /* [m] effective reflection height */
  double beta = 0.63; /* sharpness factor */

  double Lp = 0.0; /* [m] center_perturbation */
  double z_dec = 0.0; /* [m] height_decrease */
  double sig_per = 0.0; /* [m] sigma_perturbation*/

  std::complex <double> *Er0 = new std::complex <double> [Nth+1 - PML_L];

  calc_fdtd(z_prime, beta, Lp, z_dec, sig_per, Er0)

  /* Er0 の振幅と位相を出力 */
}
```

## 入江君は以下に記入して下さい。やったことなど。
