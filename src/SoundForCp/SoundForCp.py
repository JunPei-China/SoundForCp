# -*- coding: utf-8 -*-

# Cp计算小程序
import os,yaml
import numpy as np
from scipy import integrate
import scipy.constants as C
import csv

# 定义全局变量
pi = C.pi
h = C.Planck
hbar = h/(2*pi)
R = C.R
k_B = C.k
N_A = C.N_A

class HeatCapacityDebye(object):
    
    def __init__(self,name):
#        self.__sound_velocity_l = logitudinal_sound_velocity
#        self.__sound_velocity_t = transverse_sound_velocity
#        self.__number_atoms = number_atoms_value
#        self.__volume_cell = volume_cell_value
        self.__sample_name = name
        self.__debye_temperature_set = 0
        
        
    
    @property
    def sample_name(self):
        return self.__sample_name
    
    # 读写纵波声速
    @property
    def sound_velocity_l(self):
        return self.__sound_velocity_l
    @sound_velocity_l.setter
    def sound_velocity_l(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("logitudinal sound velocity must be a number")
        self.__sound_velocity_l = value

    # 读写横波声速
    @property
    def sound_velocity_t(self):
        return self.__sound_velocity_t
    @sound_velocity_t.setter
    def sound_velocity_t(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("transverse sound velocity must be a float or int")
        self.__sound_velocity_t = value
                
    #读写单胞中原子数目
    @property
    def number_atoms(self):
        return self.__number_atoms
    @number_atoms.setter
    def number_atoms(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("number atoms must be a float or int")
        self.__number_atoms = value        
    
    
    #读写物质的单胞体积
    @property
    def volume_cell(self):
        return self.__volume_cell
    @volume_cell.setter
    def volume_cell(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("volume_cell must be a float or int")
        self.__volume_cell = value
    
    #读写物质的温度，temperature,T均表示温度
    @property
    def temperature(self):
        return self.__temperature
    @temperature.setter
    def temperature(self,value):
#        if not (isinstance(value,float) or isinstance(value,int)):
#            raise ValueError("temperature must be a float or int")
        self.__temperature = value
    @property
    def T(self):
        return self.__temperature
    
    #读取平均声速
    @property
    def average_sound_velocity(self):
        value = np.power(1/3*(1/np.power(self.sound_velocity_l,3) \
            +2/np.power(self.sound_velocity_t,3)),-1/3)
        return value
    
    #读取德拜温度
    @property    
    def debye_temperature(self):
        if self.__debye_temperature_set == 0:
            value = hbar/k_B*np.power(6*pi*pi*self.number_atoms \
                /(self.volume_cell*1E-30),1/3)*self.average_sound_velocity
            return value
        elif self.__debye_temperature_set == 1:
            return self.__debye_temperature
    @debye_temperature.setter
    def debye_temperature(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("debye temperature must be a float or int") 
        else:
            self.__debye_temperature_set = 1
            self.__debye_temperature = value




    #振动积分-防溢出 
    def vibration_integral(self):
        value,err=integrate.quad(lambda x: np.power(x,4) \
            *np.exp(-x)/np.power((1-np.exp(-x)),2), \
            0,self.debye_temperature/self.temperature)
        return value    
    
#    #振动积分-原始公式
#    def vibration_integral(self):
#        value,err=integrate.quad(lambda x: np.power(x,4) \
#            *np.exp(x)/np.power((np.exp(x)-1),2), \
#            0,self.debye_temperature/self.temperature)
#        return value
    
    # 读写物质中的相对原子质量。
    @property
    def relative_atomic_mass(self):
        return self.__relative_atomic_mass
    @relative_atomic_mass.setter
    def relative_atomic_mass(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("relative atomic mass must be a number")
        self.__relative_atomic_mass = value
    
    # 读写物质中的原子的摩尔量，如1mol PbTe中含有2mol的Pb、Te原子。
    @property
    def atomic_mole_quantity(self):
        return self.__atomic_mole_quantity
    @atomic_mole_quantity.setter
    def atomic_mole_quantity(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("relative atomic mass must be a float or int")    
        self.__atomic_mole_quantity = value
    
    #读取1mol物质的热容，单位J/mol/K
    @property
    def heat_capacity_mol(self):
        value = 9*N_A*k_B*np.power(self.temperature/self.debye_temperature,3) \
            *self.vibration_integral()
        return value
    
    #读取1g物质的热容，单位J/g/K
    @property
    def heat_capacity_mass(self):
        return self.heat_capacity_mol*self.atomic_mole_quantity \
            /self.relative_atomic_mass

class HeatCapacityExpand(HeatCapacityDebye):
    def __init__(self,name):
        HeatCapacityDebye.__init__(self,name)
        self.__sample_name = name
        self.__adiabatic_bulk_modulus_set = 0
        self.__linear_expansion_coefficient_set =0
        self.__poisson_ratio_set = 0
        self.__gruneisen_constant_set = 0
        
    
    @property
    def sample_name(self):
        return self.__sample_name
    
    #读写样品的密度
    @property
    def density(self):
        return self.__density
    @density.setter
    def density(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("density must be a float or int")  
        self.__density = value
    
    # 计算泊松比
    @property
    def poisson_ratio(self):
        if self.__poisson_ratio_set == 0:
            vale_a = self.sound_velocity_t/self.sound_velocity_l
            value = (1-2*np.power(vale_a,2))/(2-2*np.power(vale_a,2))
            return value
        elif self.__poisson_ratio_set == 1:
            return self.__poisson_ratio
    @poisson_ratio.setter
    def poisson_ratio(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("poisson ratio must be a float or int") 
        else:
            self.__poisson_ratio_set = 1
            self.__poisson_ratio = value
    
    # 计算格林艾森常数
    @property
    def gruneisen_constant(self):
        if self.__gruneisen_constant_set == 0:
            value = (3*(1+self.poisson_ratio))/(2*(2-3*self.poisson_ratio))
            return value
        elif self.__gruneisen_constant_set == 1:
            return self.__gruneisen_constant
    @gruneisen_constant.setter
    def gruneisen_constant(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("gruneisen constant must be a float or int") 
        else:
            self.__gruneisen_constant_set = 1
            self.__gruneisen_constant = value
    
    #读取1mol物质的热容，单位J/mol/K
    @property
    def heat_capacity_mol_debye(self):
        value = 9*N_A*k_B*np.power(self.temperature/self.debye_temperature,3) \
            *self.vibration_integral()
        return value
    
    #读取1g物质的热容，单位J/g/K
    @property
    def heat_capacity_mass_debye(self):
        return self.heat_capacity_mol_debye* \
            self.atomic_mole_quantity/self.relative_atomic_mass
    
    # 读写绝热块体模量B_a
    @property
    def adiabatic_bulk_modulus(self):
        if self.__adiabatic_bulk_modulus_set == 0:
            value = (np.power(self.sound_velocity_l,2)-4/3  \
                *np.power(self.sound_velocity_t,2))*self.density/1E6
            return value
        elif self.__adiabatic_bulk_modulus_set == 1:
            return self.__adiabatic_bulk_modulus
    @adiabatic_bulk_modulus.setter
    def adiabatic_bulk_modulus(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("adiabatic bulk modulus must be a float or int") 
        else:
            self.__adiabatic_bulk_modulus_set = 1
            self.__adiabatic_bulk_modulus = value
    
    
    #计算剪切模量
    @property
    def shear_modulus(self):
        value=np.power(self.sound_velocity_t,2)*self.density/1E6
        return value
    
    #读取杨氏模量(该公式仅限于各向同性的材料)
    @property
    def Young_modulus(self):
        value= self.density*np.power(self.average_sound_velocity,2)*(3  \
            *np.power(self.sound_velocity_l,2)  \
            -4*np.power(self.sound_velocity_t,2))  \
            /(np.power(self.sound_velocity_l,2)   \
            -np.power(self.sound_velocity_t,2)) /1E6
        return value
    
    #读取线膨胀系数
    @property
    def linear_expansion_coefficient(self):
        if self.__linear_expansion_coefficient_set == 0:
            value = 1/3*self.gruneisen_constant*self.heat_capacity_mol_debye  \
                /(self.volume_cell*self.adiabatic_bulk_modulus)/N_A   \
                *self.number_atoms*1E21
            return value
        elif self.__linear_expansion_coefficient_set == 1:
            return self.__linear_expansion_coefficient
    @linear_expansion_coefficient.setter
    def linear_expansion_coefficient(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("density must be a float or int")  
        else:
            self.__linear_expansion_coefficient_set = 1
            self.__linear_expansion_coefficient = value
    
    
    #读取等温块体模量
    @property
    def isothermal_bulk_modulus(self):
        value = self.adiabatic_bulk_modulus  \
            /(1+self.linear_expansion_coefficient  \
            *self.gruneisen_constant*self.temperature)
        return value
    
    #读取考虑热膨胀的热容 J/g/K
    @property
    def heat_capacity_mass_expansion(self):
        value = 9*self.isothermal_bulk_modulus  \
            *np.power(self.linear_expansion_coefficient,2)  \
            *self.temperature/self.density*1E3
        return value

    #读取考虑热膨胀的热容 J/mol/K
    @property
    def heat_capacity_mol_expansion(self):
        return self.heat_capacity_mass_expansion \
                *self.relative_atomic_mass \
                /self.atomic_mole_quantity
    
    
    #读取热容 J/mol/K   
    @property
    def heat_capacity_mol(self):
        return self.heat_capacity_mol_debye+self.heat_capacity_mol_expansion
    
    #读取热容 J/mol/K  
    @property
    def heat_capacity_mass(self):
        return self.heat_capacity_mass_debye+self.heat_capacity_mass_expansion
      



def read():
    CurrentPath=os.getcwd()
    YamlFile=os.path.join(CurrentPath,"input.yaml")

    with open(YamlFile,"r",encoding="utf-8") as f:
        value = yaml.load(f,Loader=yaml.FullLoader)
    return value
def calculate():
    parameter = read()
    s = HeatCapacityExpand(parameter["Sample_Name"])
    s.sound_velocity_l = float(parameter["Longitudinal_Sound_Velocity"])
    s.sound_velocity_t = float(parameter["Transverse_Sound_Velocity"])
    s.density = float(parameter["Sample_Density"])
    s.number_atoms = float(parameter["Number_Atoms"])
    s.volume_cell = float(parameter["Volume_Cell"])
    start_temperature = float(parameter["Temperature"]["Start_Temperature"])
    end_temperature = float(parameter["Temperature"]["End_Temperature"])
    interval_temperature = float(parameter["Temperature"]["Interval_Temperature"])
    s.relative_atomic_mass = float(parameter["Relative_Atomic_Mass"])
    s.atomic_mole_quantity = float(parameter["Atomic_Mole_Quantity"])
    
    if "Debye_Temperature" in parameter.keys():
        s.debye_temperature = float(parameter["Debye_Temperature"])
    if "Poisson_Ratio" in parameter.keys():
        s.poisson_ratio = float(parameter["Poisson_Ratio"]) 
    if "Gruneisen_Constant" in parameter.keys():
        s.gruneisen_constant = float(parameter["Gruneisen_Constant"])         

    if "Adiabatic_Bulk_Modulus" in parameter.keys():
        s.adiabatic_bulk_modulus = float(parameter["Adiabatic_Bulk_Modulus"])
    if "Linear_Expansion_Coefficient" in parameter.keys():
        s.linear_expansion_coefficient = float(
            parameter["Linear_Expansion_Coefficient"])
    
    print("-"*36+"基本信息"+"-"*36)
    print(" "*4+"本程序由13skeleton编写,如有任何问题，请直接联系邮箱。(J.Pei@foxmail.com)")
    print("""
      """)
      
    
    print("-"*36+"输入参数"+"-"*36)
    print("样品名称",s.sample_name)
    print("纵波声速",s.sound_velocity_l)
    print("横波声速",s.sound_velocity_t)
    print("样品密度",s.density)
    print("单胞中原子数目",s.number_atoms)
    print("单胞体积",s.volume_cell)
    print("起始温度",start_temperature)
    print("终止温度",end_temperature)
    print("间隔温度",interval_temperature)
    print("相对原子质量",s.relative_atomic_mass)
    print("总的原子摩尔量",s.atomic_mole_quantity)
    if "Adiabatic_Bulk_Modulus" in parameter.keys():
        print("绝热块体模量",s.adiabatic_bulk_modulus)
    if "Linear_Expansion_Coefficient" in parameter.keys():
        print("线膨胀系数",s.linear_expansion_coefficient)
    print(" ")
    print(" ")
    print(" ")
    

    
    results_temperature = []
    for i in np.arange(start_temperature,end_temperature,interval_temperature):
        s.temperature = i
        list_for_temperature = (s.temperature,
            s.linear_expansion_coefficient,s.isothermal_bulk_modulus,
            s.heat_capacity_mol_debye,s.heat_capacity_mass_debye,
            s.heat_capacity_mol_expansion,s.heat_capacity_mass_expansion,
            s.heat_capacity_mol,s.heat_capacity_mass)
        results_temperature.append(list_for_temperature)

    with open("out.csv","w",newline="") as csvfile:
        myinput = csv.writer(csvfile)
        myinput.writerow(["输入参数"])
        myinput.writerow(["样品名称",s.sample_name])
        myinput.writerow(["纵波声速",s.sound_velocity_l])
        myinput.writerow(["横波声速",s.sound_velocity_t])
        myinput.writerow(["样品密度",s.density])
        myinput.writerow(["单胞中原子数目",s.number_atoms])
        myinput.writerow(["单胞体积",s.volume_cell])
        myinput.writerow(["起始温度",start_temperature])
        myinput.writerow(["终止温度",end_temperature])
        myinput.writerow(["间隔温度",interval_temperature])
        myinput.writerow(["相对原子质量",s.relative_atomic_mass])
        myinput.writerow(["总的原子摩尔量",s.atomic_mole_quantity])
        if "Adiabatic_Bulk_Modulus" in parameter.keys():
            myinput.writerow(["绝热块体模量",s.adiabatic_bulk_modulus])
        if "Linear_Expansion_Coefficient" in parameter.keys():
            myinput.writerow(["线膨胀系数",s.linear_expansion_coefficient])
        myinput.writerow([" "," "])
        

        
    with open("out.csv","a",newline="") as csvfile:
        myoutput = csv.writer(csvfile)
        myoutput.writerow(["#输出结果"])
        myoutput.writerow(["平均声速",s.average_sound_velocity])
        myoutput.writerow(["德拜温度",s.debye_temperature])
        myoutput.writerow(["泊松比",s.poisson_ratio])
        myoutput.writerow(["格林艾森常数",s.gruneisen_constant])
        myoutput.writerow(["剪切模量",s.shear_modulus])
        myoutput.writerow(["杨氏模量",s.Young_modulus])
        if "Adiabatic_Bulk_Modulus" not in parameter.keys():
            myoutput.writerow(["绝热块体模量",s.adiabatic_bulk_modulus])
        myoutput.writerow([" "," "])
        myoutput.writerow(["温度","线膨胀系数","等温块体模量","热容德拜项(mol)",
            "热容德拜项(质量)","热容膨胀项(mol)","热容膨胀项(质量)",
            "总热容(mol)","总热容(质量)"])
        myoutput.writerows(results_temperature)
    
    print("-"*36+"#输出结果"+"-"*36)
    print("平均声速",format(s.average_sound_velocity,".4f"))
    print("德拜温度",format(s.debye_temperature,".4f"))
    print("泊松比",format(s.poisson_ratio,".4f"))
    print("格林艾森常数",format(s.gruneisen_constant,".4f"))
    print("剪切模量",format(s.shear_modulus,".4f"))
    print("杨氏模量(GPa)",format(s.Young_modulus,".4f"))
    if "Adiabatic_Bulk_Modulus" not in parameter.keys():
        print("绝热块体模量(GPa)",format(s.adiabatic_bulk_modulus,".4f"))
    print("...")
    print("...")
    print("...")
            
    print("计算完成，输出结果请查看out.csv文件")

def SoundForCp():
    calculate()
    a = input("按任意键退出")    

if __name__ == "__main__":
    calculate()
    a = input("按任意键退出")
