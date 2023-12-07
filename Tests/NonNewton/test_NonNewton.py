import matplotlib.pyplot as plt
import numpy as np

class NonNewton():
    def __init__(self, viscosity0=1000.0, viscosity_inf=1.0, consistency_index=0.5, power_index=0.667, muC=10.0, yieldStress=180.0, criticalStrainRate=0.5):
        self.viscosity0 = viscosity0
        self.viscosity_inf = viscosity_inf
        self.consistency_index = consistency_index
        self.power_index = power_index
        self.muC = muC
        self.yieldStress = yieldStress
        self.criticalStrainRate = criticalStrainRate
    
    def computeViscosityNewtonian(self):
        return self.viscosity0
    
    def computeViscosityPowerLaw(self, strainRateNorm):
        return self.consistency_index * np.power(strainRateNorm, self.power_index - 1)
    
    def computeViscosityCrossModel(self, strainRateNorm):
        return self.viscosity_inf +  (self.viscosity0 - self.viscosity_inf) / (1 +  np.power(self.consistency_index * strainRateNorm, self.power_index))
    
    def computeViscosityCassonModel(self, strainRateNorm):
        return np.power(np.sqrt(self.muC) +  np.sqrt(self.yieldStress / strainRateNorm), 2)
    
    def computeViscosityCarreau(self, strainRateNorm):
        return self.viscosity_inf +  (self.viscosity0 - self.viscosity_inf) / (1.0 +  np.power(self.consistency_index * strainRateNorm * strainRateNorm, (1.0 - self.power_index)/2.0))
    
    def computeViscosityBingham(self, strainRateNorm):
        if strainRateNorm < self.criticalStrainRate:
            return self.viscosity0
        else:
            tau0 = self.criticalStrainRate * (self.viscosity0 - self.viscosity_inf)
            return tau0 / strainRateNorm + self.viscosity_inf
    
    def computeViscosityHerschelBulkley(self, strainRateNorm):
        if strainRateNorm < self.criticalStrainRate:
            return self.viscosity0
        else:
            # tau0 = self.criticalStrainRate * (self.viscosity0 - self.viscosity_inf)
            tau0 = self.viscosity0 * self.criticalStrainRate - self.consistency_index * pow(self.criticalStrainRate, self.power_index)
            return tau0 / strainRateNorm + self.consistency_index * pow(strainRateNorm, self.power_index - 1)


def test_NonNewton(model_name:str, viscosity0=1000.0, viscosity_inf=1.0, consistency_index=1000.0, power_index=0.667, muC=10.0, yieldStress=200.0, criticalStrainRate=20.0)->list:
    # TODO:更改参数！

    strainRateNorms = np.array([i*.1 for i in range(1,1000)])

    nonNewton = NonNewton(viscosity0, viscosity_inf, consistency_index, power_index, muC, yieldStress, criticalStrainRate)

    viscosity = []
    for strainRateNorm in strainRateNorms:
        if model_name == 'Newtonian':
            viscosity.append(nonNewton.computeViscosityNewtonian())
        elif model_name == 'PowerLaw':
            viscosity.append(nonNewton.computeViscosityPowerLaw(strainRateNorm))
        elif model_name == 'Cross':
            viscosity.append(nonNewton.computeViscosityCrossModel(strainRateNorm))
        elif model_name == 'Casson':
            viscosity.append(nonNewton.computeViscosityCassonModel(strainRateNorm))
        elif model_name == 'Carreau':
            viscosity.append(nonNewton.computeViscosityCarreau(strainRateNorm))
        elif model_name == 'Bingham':
            viscosity.append(nonNewton.computeViscosityBingham(strainRateNorm))
        elif model_name == 'HerschelBulkley':
            viscosity.append(nonNewton.computeViscosityHerschelBulkley(strainRateNorm))
        elif model_name == 'All':
            viscosity.append([
                nonNewton.computeViscosityNewtonian(),
                nonNewton.computeViscosityPowerLaw(strainRateNorm),
                nonNewton.computeViscosityCrossModel(strainRateNorm),
                nonNewton.computeViscosityCassonModel(strainRateNorm),
                nonNewton.computeViscosityCarreau(strainRateNorm),
                nonNewton.computeViscosityBingham(strainRateNorm),
                nonNewton.computeViscosityHerschelBulkley(strainRateNorm)])

    # 转置一下，方便画图
    if model_name == 'All':
        v_ = np.array(viscosity)
        viscosity = v_.transpose().tolist()
        
    return strainRateNorms, viscosity



def draw(x, y, color='r', label='', linestye='solid'):
    font = {'family' : 'Arial',
        'size'   : 20}

    # matplotlib.rc('font', **font)
    
    # ax = plt.subplot(111)
    # 设置刻度字体大小
    # plt.xticks(fontsize=20)
    # plt.yticks(fontsize=20)
    # 设置坐标标签字体大小
    # ax.set_xlabel(..., fontsize=20)
    # ax.set_ylabel(..., fontsize=20)
    # 设置图例字体大小
    # ax.legend(..., fontsize=20)

    SMALL_SIZE = 16
    MEDIUM_SIZE = 20
    BIGGER_SIZE = 22

    plt.rcParams['font.weight'] = 'bold'
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=BIGGER_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=BIGGER_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plt.xlabel('Shear rate(1/s)')
    plt.ylabel('Viscosity(Pa.s)')

    handle = plt.plot(x, y, color, label=label,linewidth=2, linestyle=linestye)
    return handle

def run_and_draw(model_name, **kwargs):
    [x,y] = test_NonNewton(model_name, **kwargs)
    plt.figure()
    [x,y] = test_NonNewton('Newtonian')
    np.savetxt(this_path + '/Newtonian.csv', np.array([x,y]).transpose())
    draw(x, y, 'y', label='Newtonian')
    plt.show()

if __name__ == "__main__":

    import os
    this_path = os.path.abspath(os.path.dirname(__file__))

    
    
    #将文件保存至文件中并且画出图
    # plt.savefig('figure.eps')


    # 测试所有模型
    plt.figure()

    [x,y] = test_NonNewton('Newtonian')
    # np.savetxt(this_path + '/Newtonian.csv', np.array([x,y]).transpose())
    draw(x, y, 'y', label='Newtonian', linestye='dashed')
    [x,y] = test_NonNewton('PowerLaw', power_index=0.667, consistency_index = 500.0) 
    # np.savetxt(this_path + '/PowerLaw.csv', np.array([x,y]).transpose())
    # draw(x, y, 'y', label='PowerLaw, power_index=0.667, consistency_index = 500.0')
    draw(x, y, 'y', label='PowerLaw(shear-thinning)')
    [x,y] = test_NonNewton('PowerLaw', power_index=1.1, consistency_index = 500.0)
    # np.savetxt(this_path + '/PowerLaw.csv', np.array([x,y]).transpose())
    # draw(x, y, 'b', label='PowerLaw, power_index=1.1, consistency_index = 500.0')
    draw(x, y, 'b', label='PowerLaw(shear-thickening)')
    [x,y] = test_NonNewton('Cross', consistency_index = 1.0, power_index = 0.667)
    # np.savetxt(this_path + '/Cross.csv', np.array([x,y]).transpose())
    # draw(x, y, 'c', label = 'Cross, consistency_index = 1.0, power_index = 0.667')
    draw(x, y, 'c', label = 'Cross')
    [x,y] = test_NonNewton('Casson', muC = 6.0, yieldStress = 90.0)
    # np.savetxt(this_path + '/Casson.csv', np.array([x,y]).transpose())
    # draw(x, y, 'g', label = 'Casson, muC = 6.0, yieldStress = 90.0')
    draw(x, y, 'g', label = 'Casson')
    [x,y] = test_NonNewton('Carreau', consistency_index = 0.1, power_index = 0.9)
    # np.savetxt(this_path + '/Carreau.csv', np.array([x,y]).transpose())
    draw(x, y, 'm', label = 'Carreau')
    [x,y] = test_NonNewton('Bingham', criticalStrainRate = 1.0)
    # np.savetxt(this_path + '/Bingham.csv', np.array([x,y]).transpose())
    draw(x, y, 'r', label = 'Bingham')
    [x,y] = test_NonNewton('HerschelBulkley', consistency_index = 60.0, power_index = 0.667, criticalStrainRate = 1.0)
    # np.savetxt(this_path + '/HerschelBulkley.csv', np.array([x,y]).transpose())
    draw(x, y, 'k', label = 'HerschelBulkley')

    # plt.title('All')
    plt.legend(loc="upper right")
    # plt.savefig(this_path +'/6NonNewtonsLineChart.png', dpi=300, bbox_inches='tight')




    # 测试单个模型
    test_single_model = False
    if test_single_model:
        # 测试Newtonian模型
        plt.figure()
        [x,y] = test_NonNewton('Newtonian')
        draw(x, y, 'y', label='Newtonian')
        plt.legend()
        # plt.savefig('./pics/demo_Newtonian.jpg', dpi=300)

        # 测试PowerLaw模型
        plt.figure()
        [x,y] = test_NonNewton('PowerLaw', power_index=0.667, consistency_index = 500.0) #通过给定不同参数来测试
        draw(x, y, 'r', label='power_index=0.667, consistency_index = 500.0')
        [x,y] = test_NonNewton('PowerLaw', power_index=1.1, consistency_index = 500.0)
        draw(x, y, 'b', label='power_index=1.1, consistency_index = 500.0')

        plt.title('PowerLaw')
        plt.legend()
        # plt.savefig('./pics/demo_PowerLaw.jpg', dpi=300)

        # 测试Cross模型
        plt.figure()
        [x,y] = test_NonNewton('Cross', consistency_index = 1.0, power_index = 0.667)
        draw(x, y, 'm', label = 'consistency_index = 1.0, power_index = 0.667')
        
        plt.title('Cross')
        plt.legend()
        # plt.savefig('./pics/demo_Cross.jpg', dpi = 300)

        # 测试Casson模型
        plt.figure()
        [x,y] = test_NonNewton('Casson', muC = 6.0, yieldStress = 90.0)
        draw(x, y, 'g', label = 'muC = 6.0, yieldStress = 90.0')

        plt.title('Casson')
        plt.legend()
        # plt.savefig('./pics/demo_Casson.jpg', dpi = 300)

        #测试Carreau模型
        plt.figure()
        [x,y] = test_NonNewton('Carreau', consistency_index = 0.1, power_index = 0.9)
        draw(x, y, 'm', label = 'consistency_index = 0.1, power_index = 0.9')

        plt.title('Carreau')    
        plt.legend()
        # plt.savefig('./pics/demo_Carreau.jpg', dpi = 300)

        #测试Bingham模型
        plt.figure()
        [x,y] = test_NonNewton('Bingham', criticalStrainRate = 1.0)
        draw(x, y, 'y', label = 'criticalStrainRate = 1.0')
        
        plt.title('Bingham')
        plt.legend()
        # plt.savefig('./pics/demo_Bingham.jpg', dpi = 300)

        #测试HerschelBulkley模型
        plt.figure()
        [x,y] = test_NonNewton('HerschelBulkley', consistency_index = 60.0, power_index = 0.667, criticalStrainRate = 1.0)
        draw(x, y, 'k', label = 'consistency_index = 60.0, power_index=0.667, criticalStrainRate = 1.0')
        plt.title('HerschelBulkley')
        plt.legend()
        # plt.savefig('./pics/demo_HerschelBulkley.jpg', dpi = 300)


    plt.show()