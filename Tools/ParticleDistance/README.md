用python脚本确定点云的粒子间距

粒子是从Houdini的points from volume采样的。按照我的经验，splish中的particleRadius应该取Houdini particle seperation的一半。否则过密会爆炸，过稀疏会很散。

这个脚本用来计算输入的点云的粒子间距是多少，以此进行判断。