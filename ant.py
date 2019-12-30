import string
import matplotlib.pyplot as plt
import numpy as np
import math

class AS():
	
	def __init__(self, m = 50, Nc = 100, a = 1, b = 6, p = 0.9, q = 50):
		self.m = m		#蚂蚁个数
		self.Nc = Nc	#迭代次数
		self.a = a		#alpha参数，信息启发因子
		self.b = b		#beta参数，期望启发因子
		self.p = p		#信息素挥发因子
		self.q = q		#信息素强度
		self.distance = np.zeros((1, 1))	#城市间距离矩阵
		self.eta = np.zeros((1, 1))		#启发函数矩阵
		self.phero = np.ones((1, 1))	#信息素浓度矩阵
		self.x = []		#城市的x坐标
		self.y = []		#城市的y坐标
		self.bestpath = []		#最佳路径
		#self.readtxt("ei151.txt")
		
	
	def readtxt(self, filename = "ei151.txt"):		#读取城市坐标
		with open(filename) as f:
			for lines in f:
				temp = lines.split()
				self.x.append(float(temp[1]))
				self.y.append(float(temp[2]))
		num = len(self.x)
		self.distance = np.zeros((num, num))		#城市距离矩阵初始化
		self.eta = np.zeros((num, num))				#启发函数矩阵初始化
		self.phero = np.ones((num, num))			#信息素浓度矩阵初始化
		for i in range(num):
			for j in range(i, num):
				self.distance[i][j] = self.distance[j][i] = math.sqrt(((self.x[i]-self.x[j])**2)+((self.y[i]-self.y[j])**2))		#欧氏距离矩阵
				if(self.distance[i][j] != 0):
					self.eta[i][j] = self.eta[j][i] = 1 / self.distance[i][j]		#启发信息矩阵

	def iterator(self):
		iternum = 0		#初始化迭代次数
		while iternum < self.Nc:
			delta = self.p * self.phero		#初始化信息素变换矩阵
			distance_ant = []		#迭代过程中的所有蚂蚁周游距离
			path_bests = []		#存储最佳路径
			for i in range(self.m):
				#start = np.random.randint(len(self.x))		#随机化蚂蚁初始位置
				start = 0
				path = [start]		#路径
				tabu = [start]		#禁忌表
				d = 0		#周游距离
				for j in range(len(self.x)-1):
					non_tabu = list(set(range(len(self.x))) - set(tabu))
					movep = np.zeros(len(self.x)-len(tabu))		#转移概率矩阵
					for k in range(len(self.x)-len(tabu)):
						movep[k] = (self.phero[path[-1]][non_tabu[k]] ** self.a) * (self.eta[path[-1]][non_tabu[k]] ** self.b) 		#转移概率
					p = movep / sum(movep)
					r = np.random.rand()		#轮盘赌算法
					'''
					while True:
						r = np.random.rand()
						index = np.where(p > r)[0]
						if len(index) > 0:
							l = non_tabu[index[0]]
							break
					'''
					accumulator = 0.0
					i = 0
					for i in range(len(p)):
						accumulator += p[i]
						if accumulator >= r:
							break
					l = non_tabu[i]

					d = d + self.distance[path[-1]][l]
					path.append(l)
					tabu.append(l)
				d = d + self.distance[path[-1]][start]		#路径形成闭环
				path.append(start)
				distance_ant.append(d)
				path_bests.append(path)
				for j in range(len(path)-1):
					delta[path[j]][path[j+1]] += self.q / d		#蚁环增量模型
					#delta[path[j]][path[j + 1]] += self.q * self.distance[path[j]][path[j + 1]]/ d		#蚁恒增量模型
					delta[path[j+1]][path[j]] = delta[path[j]][path[j+1]]
			print(min(distance_ant))
			index = distance_ant.index(min(distance_ant))
			print(path_bests[index])
			for j in range(len(path_bests[index]) - 1):
				delta[path_bests[index][j]][path_bests[index][j + 1]] += self.q  / min(distance_ant)		#强化最优路径信息素
				#delta[path_bests[index][j+1]][path_bests[index][j]] = delta[path_bests[index][j]][path_bests[index][j + 1]]
			self.phero = delta		#更新信息素矩阵
			self.bestpath = path_bests[index]		#更新最优路径
			iternum += 1



	def display(self):		#可视化路径
		'''
		x = [[1, 3], [2, 5]]
		y = [[4, 7], [6, 3]]
		for i in range(len(x)):
			plt.plot(x[i], y[i], color='r')
			plt.scatter(x[i], y[i], color='b')
		plt.show()
		plt.scatter(self.x, self.y, s=10)
		'''
		path = []
		for j in range(len(self.bestpath)):
			path.append([self.x[self.bestpath[j]],self.y[self.bestpath[j]]])

		for j in range(len(self.bestpath)-1):
			plt.scatter(path[j][0], path[j][1], s=10, color = 'b')
			plt.plot([path[j][0], path[j+1][0]],[path[j][1], path[j+1][1]], color = 'g')
		print(path)
		plt.show()


if __name__ == '__main__':
	asc = AS()
	asc.readtxt()
	asc.iterator()
	asc.display()
