import numpy as np
from PIL import Image
from sklearn.model_selection import KFold
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.model_selection import GridSearchCV
from sklearn.svm import SVC
from sklearn.ensemble import RandomForestClassifier
import math
import time
import os, sys





def uniform(src):
    table = [0] * 256
    temp = 1
    for i in range(256):
        count = 0
        t = '{:08b}'.format(i)
        for x in range(8):
            if t[x] != t[(x + 1) % 8]:
                count += 1
        if count < 3:
            table[i] = temp
            temp += 1
    #print(table)
    src = np.array(src)
    height = src.shape[0]
    width = src.shape[1]
    dst = src.copy()
    for x in range(1, width - 1):
        for y in range(1, height - 1):
            dst[y,x] = table[src[y,x]]
    dst = Image.fromarray(dst)
    return dst



def LBP(src):
    src = np.array(src)
    height = src.shape[0]
    width = src.shape[1]
    dst = src.copy()

    lbp_value = np.zeros((1, 8), dtype=np.uint8)
    neighbours = np.zeros((1, 8), dtype=np.uint8)

    for x in range(1, width - 1):
        for y in range(1, height - 1):
            neighbours[0, 0] = src[y - 1, x - 1]
            neighbours[0, 1] = src[y - 1, x]
            neighbours[0, 2] = src[y - 1, x + 1]
            neighbours[0, 3] = src[y, x - 1]
            neighbours[0, 4] = src[y, x + 1]
            neighbours[0, 5] = src[y + 1, x - 1]
            neighbours[0, 6] = src[y + 1, x]
            neighbours[0, 7] = src[y + 1, x + 1]

            center = src[y, x]

            for i in range(8):
                if neighbours[0, i] > center:
                    lbp_value[0, i] = 1
                else:
                    lbp_value[0, i] = 0

            lbp = lbp_value[0, 0] * 1 + lbp_value[0, 1] * 2 + lbp_value[0, 2] * 4 + lbp_value[0, 3] * 8 \
                  + lbp_value[0, 4] * 16 + lbp_value[0, 5] * 32 + lbp_value[0, 6] * 64 + lbp_value[0, 0] * 128

            dst[y, x] = lbp
    dst = Image.fromarray(dst)
    return dst


def circular_LBP(src, radius, n_points):
    src = np.array(src)
    height = src.shape[0]
    width = src.shape[1]
    dst = src.copy()
    src.astype(dtype=np.float32)
    dst.astype(dtype=np.float32)

    neighbours = np.zeros((1, n_points), dtype=np.uint8)
    lbp_value = np.zeros((1, n_points), dtype=np.uint8)
    for x in range(radius, width - radius - 1):
        for y in range(radius, height - radius - 1):
            lbp = 0.
            # 先计算共n_points个点对应的像素值，使用双线性插值法
            for n in range(n_points):
                theta = float(2 * np.pi * n) / n_points
                x_n = x + radius * np.cos(theta)
                y_n = y - radius * np.sin(theta)

                # 向下取整
                x1 = int(math.floor(x_n))
                y1 = int(math.floor(y_n))
                # 向上取整
                x2 = int(math.ceil(x_n))
                y2 = int(math.ceil(y_n))

                # 将坐标映射到0-1之间
                tx = np.abs(x - x1)
                ty = np.abs(y - y1)

                # 根据0-1之间的x，y的权重计算公式计算权重
                w1 = (1 - tx) * (1 - ty)
                w2 = tx * (1 - ty)
                w3 = (1 - tx) * ty
                w4 = tx * ty

                # 根据双线性插值公式计算第k个采样点的灰度值
                neighbour = src[y1, x1] * w1 + src[y2, x1] * w2 + src[y1, x2] * w3 + src[y2, x2] * w4

                neighbours[0, n] = neighbour

            center = src[y, x]

            for n in range(n_points):
                if neighbours[0, n] > center:
                    lbp_value[0, n] = 1
                else:
                    lbp_value[0, n] = 0

            for n in range(n_points):
                lbp += lbp_value[0, n] * 2 ** n

            # 转换到0-255的灰度空间，比如n_points=16位时结果会超出这个范围，对该结果归一化
            dst[y, x] = int(lbp / (2 ** n_points - 1) * 255)
    dst = Image.fromarray(dst)
    return dst


def value_rotation(num):
    value_list = np.zeros((8), np.uint8)
    temp = int(num)
    value_list[0] = temp
    for i in range(7):
        temp = ((temp << 1) | int(temp / 128)) % 256
        value_list[i + 1] = temp
    # print value_list
    return np.min(value_list)


def rotation_invariant_LBP(src):
    src = np.array(src)
    height = src.shape[0]
    width = src.shape[1]
    dst = src.copy()

    lbp_value = np.zeros((1, 8), dtype=np.uint8)
    neighbours = np.zeros((1, 8), dtype=np.uint8)
    for x in range(1, width - 1):
        for y in range(1, height - 1):
            neighbours[0, 0] = src[y - 1, x - 1]
            neighbours[0, 1] = src[y - 1, x]
            neighbours[0, 2] = src[y - 1, x + 1]
            neighbours[0, 3] = src[y, x - 1]
            neighbours[0, 4] = src[y, x + 1]
            neighbours[0, 5] = src[y + 1, x - 1]
            neighbours[0, 6] = src[y + 1, x]
            neighbours[0, 7] = src[y + 1, x + 1]

            center = src[y, x]

            for i in range(8):
                if neighbours[0, i] > center:
                    lbp_value[0, i] = 1
                else:
                    lbp_value[0, i] = 0

            lbp = lbp_value[0, 0] * 1 + lbp_value[0, 1] * 2 + lbp_value[0, 2] * 4 + lbp_value[0, 3] * 8 \
                  + lbp_value[0, 4] * 16 + lbp_value[0, 5] * 32 + lbp_value[0, 6] * 64 + lbp_value[0, 0] * 128

            # 旋转不变值
            dst[y, x] = value_rotation(lbp)
    dst = Image.fromarray(dst)
    return dst


def pre_data():
    for person in range(1, 41):
        for face in range(1, 11):
            #I = Image.open('D:\\pyproject\\venv\\ORL\\ORL\\s' + str(person) + '_' + str(face) + '.bmp').convert('L')
            I = Image.open('D:\\pyproject\\venv\\ORL\\CLBP\\s' + str(person) + '_' + str(face) + '.bmp').convert('L')
            #I.save('D:\\pyproject\\venv\\ORL\\GRAY\\s' + str(person) + '_' + str(face) + '.bmp')
            # I.save('11.PNG')
            #I_lbp = LBP(I)
            #I_lbp.save('D:\\pyproject\\venv\\ORL\\LBP\\s' + str(person) + '_' + str(face) + '.bmp')
            #I_lbp_circle = circular_LBP(I, 2, 16)
            # I_lbp_circle.show()
            #I_lbp_circle.save('D:\\pyproject\\venv\\ORL\\CLBP\\s' + str(person) + '_' + str(face) + '.bmp')
            #I_lbp_rotate_invariant = rotation_invariant_LBP(I)
            # I_lbp_rotate_invariant.show()
            #I_lbp_rotate_invariant.save('D:\\pyproject\\venv\\ORL\\RILBP\\s' + str(person) + '_' + str(face) + '.bmp')
            I_uniform = uniform(I)
            I_uniform.save('D:\\pyproject\\venv\\ORL\\ULBP\\s' + str(person) + '_' + str(face) + '.bmp')


def get_picture():
    for person in range(1, 41):
        for face in range(1, 11):
            #I = Image.open('D:\\pyproject\\venv\\ORL\\ORL\\s' + str(person) + '_' + str(face) + '.bmp').convert('L')
            #I = Image.open('D:\\pyproject\\venv\\ORL\\RILBP\\s' + str(person) + '_' + str(face) + '.bmp').convert('L')
            #I = Image.open('D:\\pyproject\\venv\\ORL\\CLBP\\s' + str(person) + '_' + str(face) + '.bmp').convert('L')
            #I = Image.open('D:\\pyproject\\venv\\ORL\\LBP\\s' + str(person) + '_' + str(face) + '.bmp').convert('L')
            I = Image.open('D:\\pyproject\\venv\\ORL\\ULBP\\s' + str(person) + '_' + str(face) + '.bmp').convert('L')
            all_data_set.append(list(I.getdata()))
            all_data_label.append(person)
    print("read text...")
    print(np.array(all_data_set))
    print(np.array(all_data_set).shape,np.array(all_data_set).dtype)

def RandomForest():
    kf = KFold(n_splits=10, shuffle=True)
    precision_average = 0.0
    clf = RandomForestClassifier(n_estimators=15)
    #clf = SVC()
    for train, test in kf.split(x):
        #print("fit " + str(train) + "...")
        clf = clf.fit(x[train], y[train])
        test_pred = clf.predict(x[test])

        precision = 0
        for i in range(0, len(y[test])):
            if y[test][i] == test_pred[i]:
                precision = precision + 1
        precision_average += float(precision) / len(y[test])
    precision_average = precision_average / 10
    print(u"准确率为" + str(precision_average))


all_data_set = []
all_data_label = []
if __name__ == '__main__':

    #pre_data()

    get_picture()
    #pca = PCA(n_components=59).fit(np.array(all_data_set))
    #all_data_pca = pca.transform(all_data_set)
    #print(np.array(all_data_pca))
    x = np.array(all_data_set)
    y = np.array(all_data_label)
    time_start = time.time()
    RandomForest()
    time_end = time.time()
    print('time cost', time_end - time_start)

