#n2_simulation_output
#シミュレーションデータからcos, 2D作成	

from scipy.special import sph_harm as SH
import numpy as np
import scipy
import os
import math
import glob

#定数
t_start = 0 #初期時間
t_step = 0 #時間刻み、シミュレーション：3.33333E-15or1.33333e-15or1.16667E-15	
angle = np.arange(0, 180, 0.5) #角度刻み
Jmax = 20 #simulation時のnonadiabatichにあるJ値
B0 = 1.989581 #分子ごとの定数, N2パターン
De = 5.76E-6
c0 = 29979250000 

#ファイル読み込み


#角度分布計算
#M=0ver.
for j in range(0,Jmax-2): #J=0
	memo[j]={}
	Erot = B0*j*(j+1)-De*j*j*(j+1)*(j+1)
	Erot2 = Erot*c0
	for m in range(-j, j+1, 2):
		memo[j][m] = np.zeros((A), dtype = np.complex)
		for ang in angle:
			temp = SH(m,j,phi,ang/(A/2.)*PI)

for j in range(1,Jmax-1): #J=1
	memo[j]={}
	Erot = B0*j*(j+1)-De*j*j*(j+1)*(j+1)
	Erot2 = Erot*c0
	for m in range(-j, j+1, 2)

#M=1,-1ver.
for j in range(1,Jmax-1): #J=1
	memo[j]={}
	Erot = B0*j*(j+1)-De*j*j*(j+1)*(j+1)
	Erot2 = Erot*c0
	for m in range(-j, j+1, 2)

#ファイル出力
#レーザーなまらせ
#加重平均行列作成
#２Dカーペット計算パートcos2計算パート
#cos2計算パート
#cos2出力パート、２Dカーペット
#polarplot出力パート
