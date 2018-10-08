#n2_simulation_output
#シミュレーションデータからcos, 2D作成	

from scipy.special import sph_harm as SH
import numpy as np
import scipy
import os
import glob

#定数
t_start = 0 #初期時間
t_step = 0 #時間刻み、シミュレーション：3.33333E-15or1.33333e-15or1.16667E-15	
angle = np.arange(0.5, 180.5, 0.5) #角度刻み
Jmax = 30 #simulation時のnonadiabatichにあるJ値
B0 = 1.989581 #分子ごとの定数, N2パターン
De = 5.76E-6
c0 = 29979250000 
phi = 0.0
PI = np.pi

#angle_count = 360/0.5 + 1
#print(angle_count)

#なまらせパート
region = 120 #前後の読み込み範囲、30で48、100で120、170で137に１を足す
point = region + 1 #初期時間、ファイル番号、region+1
delta_t = 1.33333e-15 #シミュレーションの刻み幅、30=3.33333E-15、100=1.33333e-15、170=1.16667E-15
large_t = 12e-15 #画像刻み幅、12or20or30fs
kizami = int(large_t / delta_t) #繰り返す刻み間隔、切り捨て整数化
#print(kizami), 3
point_time = point * delta_t * 1e+15 #時間計算(初期)、フェムト単衣


#length_memo = length(Jstart, Jmax, Mstart)
def length(Jstart, Jmax, Mstart):
	list1 = []
	for j in range(Jstart, Jmax, 2):
		for m in range(-(j-Mstart), j+1, 2):
			list1.append(0)
	length_memo = len(list1)*2
	return(length_memo)


def sph_harm(Jstart, Jmax, Mstart, B0, De, phi, PI, c0, memo): #球面調和関数計算
	for j in range(Jstart, Jmax, 2):
		Erot = B0*j*(j+1)-De*j*j*(j+1)*(j+1)
		for m in range(-(j-Mstart), j+1, 2):
			Erot_list.append(Erot*c0)
			Erot_list.append(Erot*c0)
			harm_list = [SH(m,j,phi,(ang/(180.))*PI) for ang in np.arange(0, 180.5, 0.5)]
			harm_list.extend(harm_list[::-1][1:-1]) #0から359.5度まで
			harm_list2 = [1j*SH(m,j,phi,(ang/(180.))*PI) for ang in np.arange(0, 180.5, 0.5)]
			harm_list2.extend(harm_list2[::-1][1:-1]) #0から359.5度まで
			memo = np.vstack((memo, np.array(harm_list)))
			memo = np.vstack((memo, np.array(harm_list2))) #球面調和関数
	return(memo, Erot_list)
			#print(len(memo))
'''
def energy(Erot_list, c0, PI, Tlist):
	energy1 = np.dot(Tlist, np.array([Erot_list])) #t行jm列の行列
	energy2 = np.exp(-1.0j*2*PI*energy1)
'''

def weight_ave(region, delta_t): #加重平均行列作成
	gauss_sigma = 80e-15 / (2 * np.sqrt(2 * np.log(2))) #シグマの値
	ave = [] #加重平均リスト
	ave_sum = 0 #加重平均要素和
	for r in range(-region, region+1):
		#print(r)
		y = (1 / np.sqrt(2 * np.pi * gauss_sigma ) ) * np.exp(-(r * delta_t) ** 2 / (2 * (gauss_sigma ** 2)) ) #ガウス関数
		ave.append(y)
		ave_sum = ave_sum + y
	ave1 = np.array(ave)
	ave2 = ave1.reshape(region*2+1, 1) #転置
	return(ave2)

#角度分布計算
kensaku = glob.glob('[0-1][0-2]terms.dat') #00~12までを検索, 00.10.11.12
for filename in kensaku:
	Jpop = np.loadtxt(open(filename), skiprows = 1) #１行飛ばす
	Tlist = Jpop[:,-1].reshape(len(Jpop[:,-1]),1)
	Jpop = Jpop[:,0:-1] #時間削除
	if '00' in filename:
		Erot_list=[]
		Jstart = 0
		Mstart = 0
		#length_memo = length(Jstart, Jmax, Mstart)
		memo = np.zeros((0, 720), dtype=np.complex)
		memo = sph_harm(Jstart, Jmax, Mstart, B0, De, phi, PI, c0, memo)[0] #球面調和関数計算
		Erot_list = sph_harm(Jstart, Jmax, Mstart, B0, De, phi, PI, c0, memo)[1] #エネルギー
		#エネルギー項の計算
		energy1 = np.dot(Tlist, np.array([Erot_list])) #t行jm列の行列
		energy2 = np.exp(-1.0j*2*PI*energy1)
		#行列計算、角度分布計算
		pop_cal = np.dot(Jpop, memo) #行列の掛け算
		#共役を掛ける
		carpet00 = pop_cal * np.conj(pop_cal)
	elif '10' in filename:
		Jpop10 = Jpop
	elif '11' in filename:
		Erot_list=[]
		Jstart = 1
		#length_memo = length(Jstart, Jmax, Mstart)
		memo = np.zeros((0, 720), dtype=np.complex)
		sph_harm(Jstart, Jmax, Mstart, B0, De, phi, PI, c0, memo) #球面調和関数計算
		#行列計算、角度分布計算
		pop_cal = np.dot(Jpop, memo) #行列の掛け算
		carpet11 = pop_cal * np.conj(pop_cal)
	elif '12' in filename:
		Jpop = Jpop10/9 + Jpop/9
		Erot_list=[]
		Jstart = 1
		Mstart = 0 #初期値はJM=1-1,11
		#length_memo = length(Jstart, Jmax, Mstart)
		memo = np.zeros((0, 720), dtype=np.complex) #角度が0.5°刻みのため、720コ
		sph_harm(Jstart, Jmax, Mstart, B0, De, phi, PI, c0, memo) #球面調和関数計算
		#行列計算、角度分布計算
		pop_cal = np.dot(Jpop, memo) #行列の掛け算
		carpet1012 = pop_cal * np.conj(pop_cal)

carpet = carpet00 * 2/3 + carpet11 / 9 + carpet1012 #carpet完成

#ファイル出力
#レーザーなまらせ
weight = weight_ave(region, delta_t)#加重平均行列作成


#２Dカーペット計算パートcos2計算パート
#cos2計算パート
#cos2出力パート、２Dカーペット
#polarplot出力パート


'''
#M=0ver.
Erot_list=[]
memo = np.zeros((0, 720))
Jstart = 0
Mstart = 0
sph_harm(Jstart, Jmax, Mstart, B0, De, phi, PI) #球面調和関数計算

exp(-1.0j*Erot*c0*Tnow)


#角度分布計算
#dic作成、球面調和関数
#M=0ver.
for j in range(0,Jmax,2): #J=0
	memo1={}
	Erot1={}
	Erot = B0*j*(j+1)-De*j*j*(j+1)*(j+1)
	Erot1[j] = Erot*c0
	for m in range(-j, j+1, 2):
		memo1[(j, m)] = np.zeros((721), dtype = np.complex)
		ang = 0
		temp = SH(m,j,phi,(ang/(180.))*PI)
		memo1[(j, m)][ang] = temp
		for ang in angle:
			ang2 = 360-ang
			temp = SH(m,j,phi,(ang/(180.))*PI)
			memo1[(j, m)][ang] = temp
			memo1[(j, m)][ang2] = temp


for j in range(1,Jmax,2): #J=1
	memo2={}
	Erot2={}
	Erot = B0*j*(j+1)-De*j*j*(j+1)*(j+1)
	Erot2[j] = Erot*c0
	for m in range(-(j-1), j+1, 2):
		memo2[(j, m)] = np.zeros((721), dtype = np.complex)
		ang = 0
		temp = SH(m,j,phi,(ang/(180.))*PI)
		memo2[(j, m)][ang] = temp
		for ang in angle:
			ang2 = 360-ang
			temp = SH(m,j,phi,(ang/(180.))*PI)
			memo2[(j, m)][ang] = temp
			memo2[(j, m)][ang2] = temp


#M=1,-1ver.
for j in range(1,Jmax, 2): #J=1
	memo3={}
	Erot3={}
	Erot = B0*j*(j+1)-De*j*j*(j+1)*(j+1)
	Erot3[j] = Erot*c0
	for m in range(-j, j+1, 2):
		memo3[(j, m)] = np.zeros((721), dtype = np.complex)
		ang = 0
		temp = SH(m,j,phi,(ang/(180.))*PI)
		memo3[(j, m)][ang] = temp
		for ang in angle:
			ang2 = 360-ang
			temp = SH(m,j,phi,(ang/(180.))*PI)
			memo3[(j, m)][ang] = temp
			memo3[(j, m)][ang2] = temp

#ファイル読み込み
kensaku = glob.glob('[0-1][0-2]terms.dat') #00~12までを検索, 00.10.11.12
gyo3 = np.vstack((gyo3, np.loadtxt(open(filename))[:-1]))
[:,-1] #時間抽出
b[:,0:-1] #時間削除

#始状態J=1
#Mが奇数のものだけ先に重ね合わせる
'''







