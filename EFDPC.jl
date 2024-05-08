# 导入相关包：
using Statistics
using MAT
using Base
using LinearAlgebra

# 加载原始图像数据
data = matread("indian_pines.mat")      # 从MAT文件中读取数据
Img = data["indian_pines"][:, :, vcat(4:102, 113:147, 166:216)]  # 选择特定波段的数据
Nr, Nc, Nb = size(Img)                  # 获取数据的维度信息
Img_matrix = reshape(Img, Nr * Nc, Nb)  # 不需要空间信息，直接将高光谱图像数据转换为二维矩阵进行处理

# 对波段进行预处理
for i = 1 : Nr * Nc
    Img_matrix[i, :] .= (Img_matrix[i, :] .- minimum(Img_matrix[i, :])) ./ (maximum(Img_matrix[i, :]) - minimum(Img_matrix[i, :]))
end

# 第一步：构建相似度矩阵
Dist_matrix = zeros(Nb, Nb)     # 初始化相似度矩阵
for i = 1 : Nb-1    # 循环遍历每个波段
    Vi = Img_matrix[:, i]
    for j = i+1 : Nb
        Vj = Img_matrix[:, j]
        Dist_matrix[i, j] = norm(Vi - Vj)       # 计算2范数
        Dist_matrix[j, i] = Dist_matrix[i, j]   # 距离矩阵是对称的
    end
end
Dist_matrix ./= Nb      # 缩放相似度矩阵获得距离矩阵

# 第二步：计算局部密度
# 计算dc
percent = 0.02    # 选择倒数2%作为dc值
position = Int(round(Nb * (Nb - 1) / 2 * percent))  # 找到dini的位置2% × L × (L − 1)
temp = Dist_matrix[tril(Dist_matrix) .!= 0]         # 从 Dist_matrix 中提取出所有下三角部分的非零元素，并存储到数组 temp 中。
sda = sort(temp)        # 从小到大排序，sort() 函数用于对数组或向量进行排序，并返回排序后的新数组或向量。
dini = sda[position]    # 找到dini
k = 10  # 需要选择的波段数
dc = dini / exp(k / Nb)  # 计算dc值

# 计算 𝜌(rho) 因子
rho = zeros(Nb)
for i = 1 : Nb
    for j = 1 : Nb
        if i != j
            rho[i] += exp(-(Dist_matrix[i, j] / dc)^2)  # 根据距离计算rho
        end
    end
end

# 第三步：计算相对距离
# 计算 𝛿(delta) 因子
ordrho = sortperm(rho, rev=true)    # 获得𝜌从大到小排序序号，sortperm() 函数用于返回数组或向量排序后的索引数组（排列的是索引），而不是直接返回排序后的值。
delta = zeros(Nb)       # 初始化delta数组
delta[ordrho[1]] = -1   # 处理 𝜌 值最大的数据点：密度最高的样本不存在比其密度更高的点，先设置为最小值

for i = 2 : Nb
    delta[ordrho[i]] = floatmax(Float64)
    for j = 1 : i-1     # 在𝜌比i大的波段中，遍历寻找一个距离最近的波段j
        delta[ordrho[i]] = min(delta[ordrho[i]], Dist_matrix[ordrho[i], ordrho[j]])
    end
end
delta[ordrho[1]] = maximum(delta)   # 生成 𝜌 值最大数据点的delta值

# 归一化因子
rho = (rho .- minimum(rho)) ./ (maximum(rho) - minimum(rho))
delta = (delta .- minimum(delta)) ./ (maximum(delta) - minimum(delta))

# 第四步：计算聚类得分
gamma = rho .* delta .* delta

# 找到选择的波段
order_band = sortperm(gamma, rev=true)
C = order_band[1:k]