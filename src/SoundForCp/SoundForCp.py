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
    
    #读写样品C12和C44的关系
    @property
    def elastic_constants_condition(self):
        return self.__elastic_constants_condition
    @elastic_constants_condition.setter
    def elastic_constants_condition(self,value):
        if not (isinstance(value,float) or isinstance(value,int)):
            raise ValueError("density must be a float or int")
        self.__elastic_constants_condition = value
    
    #读写样品C12和C44的关系
    @property
    def sound_velocity_rms(self):
        value = np.power((np.power(self.sound_velocity_l,2) \
                + 2*np.power(self.sound_velocity_t,2))/3,1/2)
        return value
        
    
    # 读写绝热块体模量B_a
    @property
    def adiabatic_bulk_modulus(self):
        if self.__adiabatic_bulk_modulus_set == 0:
            if self.elastic_constants_condition == 0:
                value = (1+self.poisson_ratio)/(2-3*self.poisson_ratio)*self.density \
                        *np.power(self.sound_velocity_rms,2) / 1E6
            elif self.elastic_constants_condition == 1:
                value = self.density *np.power(self.sound_velocity_rms,2) / 1E6
            elif self.elastic_constants_condition == 2:
                value = (np.power(self.sound_velocity_l,2)-4/3  \
                    *np.power(self.sound_velocity_t,2))*self.density/1E6
            else:
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
            /(1+3*self.linear_expansion_coefficient  \
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
    s.elastic_constants_condition = float(parameter["Elastic_Constants_Condition"])
    
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
    
    print("+"*37+"  SoundForCp  "+"+"*37)
    print("+"+" "*35+"version: 0.04"+" "*38+"+")
    print("+"+" "*27+"Developed by Jing-Feng Li's Group"+" "*26+"+")
    print("+"*88)  
    print(" ")
    print("-"*34+"  Input Parameters  "+"-"*34)
    print("Sample Name",s.sample_name)
    print("LongitudinalSound Velocity",s.sound_velocity_l)
    print("Transverse Sound Velocity",s.sound_velocity_t)
    print("Sample Density",s.density)
    print("Number Atoms in Unit Cell",s.number_atoms)
    print("Unit Cell Volume",s.volume_cell)
    print("Start Temperature",start_temperature)
    print("End Temperature",end_temperature)
    print("Interval Temperature",interval_temperature)
    print("Relative Atomic Mass",s.relative_atomic_mass)
    print("Atomic Mole Quantity",s.atomic_mole_quantity)
    if s.elastic_constants_condition == 0:
        print("C12 ≠ C44")
    elif s.elastic_constants_condition == 1:
        print("C12 = C44")
    else:
        print("polycrystals without prefered orientation.")
        
    if "Adiabatic_Bulk_Modulus" in parameter.keys():
        print("Adiabatic Bulk Modulus",s.adiabatic_bulk_modulus)
    if "Linear_Expansion_Coefficient" in parameter.keys():
        print("Linear Expansion Coefficient",s.linear_expansion_coefficient)
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
        myinput.writerow(["Input Parameters"])
        myinput.writerow(["Sample Name",s.sample_name])
        myinput.writerow(["Longitudinal Sound Velocity",s.sound_velocity_l])
        myinput.writerow(["Transverse Sound Velocity",s.sound_velocity_t])
        myinput.writerow(["Sample Density",s.density])
        myinput.writerow(["Number Atoms in Unit Cell",s.number_atoms])
        myinput.writerow(["Volume Unit Cell",s.volume_cell])
        myinput.writerow(["Start Temperature",start_temperature])
        myinput.writerow(["End Temperature",end_temperature])
        myinput.writerow(["Interval Temperature",interval_temperature])
        myinput.writerow(["Relative Atomic Mass",s.relative_atomic_mass])
        myinput.writerow(["Atomic Mole Quantity",s.atomic_mole_quantity])
        myinput.writerow(["Elastic Constants Condition",s.elastic_constants_condition])
        if "Adiabatic_Bulk_Modulus" in parameter.keys():
            myinput.writerow(["Adiabatic Bulk Modulus",s.adiabatic_bulk_modulus])
        if "Linear_Expansion_Coefficient" in parameter.keys():
            myinput.writerow(["Linear Expansion Coefficient",s.linear_expansion_coefficient])
        myinput.writerow([" "," "])
                
    with open("out.csv","a",newline="") as csvfile:
        myoutput = csv.writer(csvfile)
        myoutput.writerow(["#Output Results"])
        myoutput.writerow(["Average Sound Velocity",s.average_sound_velocity])
        myoutput.writerow(["Debye Temperature",s.debye_temperature])
        myoutput.writerow(["Poisson Ratio",s.poisson_ratio])
        myoutput.writerow(["Gruneisen Constant",s.gruneisen_constant])
        myoutput.writerow(["Shear Modulus",s.shear_modulus])
        myoutput.writerow(["Young's Modulus",s.Young_modulus])
        if "Adiabatic_Bulk_Modulus" not in parameter.keys():
            myoutput.writerow(["Adiabatic Bulk Modulus",s.adiabatic_bulk_modulus])
        myoutput.writerow([" "," "])
        myoutput.writerow(["Temperature(K)","Linear_Expansion_Coefficient(/K)","Isothermal Bulk Modulus","Isobaric Heat Capacity of Debye Term(J/mol/K)",
            "Isobaric Heat Capacity of Debye Term(J/g/K)","Isobaric Heat Capacity of Expansion Term(J/mol/K)","Isobaric Heat Capacity of Expansion Term(J/g/K)",
            "Isobaric Heat Capacity(J/mol/K)","Isobaric Heat Capacity(J/g/K)"])
        myoutput.writerows(results_temperature)
    
    print("-"*35+"  Output Results  "+"-"*35)
    print("Average Sound Velocity(m/s)",format(s.average_sound_velocity,".4f"))
    print("Debye Temperature(K)",format(s.debye_temperature,".4f"))
    print("Poisson Ratio",format(s.poisson_ratio,".4f"))
    print("Gruneisen Constant",format(s.gruneisen_constant,".4f"))
    print("Shear Modulus (GPa)",format(s.shear_modulus,".4f"))
    print("Young Modulus (GPa)",format(s.Young_modulus,".4f"))
    if "Adiabatic_Bulk_Modulus" not in parameter.keys():
        print("Adiabatic Bulk Modulus (GPa)",format(s.adiabatic_bulk_modulus,".4f"))
    print("...")
    print("...")
    print("...")
            
    print("Calculation was completed，Details can be seen in out.csv file.")
    
    print(" ")
    print("+"*40+"  Tips  "+"+"*40)
    print("+"+" "*4+"This program was developed by Jun Pei & Jing-Feng Li from Tsinghua University. If "+"+")
    print("+"+"you have any question, please feel free to contact us."+" "*32+"+")
    print("+"+"J.Pei@foxmail.com & Jingfeng@mail.tsinghua.edu.cn."+" "*36+"+")
    print("+"+" "*86+"+")
    print("+"+"Reference:"+" "*76+"+")
    print("+"+" "*4+"1. Pei, J., Li, H., Zhuang, H.-L., Dong, J., Cai, B., Hu, H., Li, J.-W., Jiang, Y."+"+")
    print("+"+" "*4+", Su, B.,Zhao, L.-D., & Li, J.-F. A sound velocity method for determining isobaric"+"+")
    print("+"+" "*4+"specific heat capacity. InfoMat,2022, e12372. https://doi.org/10.1002/inf2.12372  "+"+")
    print("+"+" "*86+"+")
    print("+"+"If you have used SoundForCp, please cite the above article."+" "*27+"+")
    print("+"*88)
    

def SoundForCp():
    calculate()
    a = input("press any keys to quit...")    

if __name__ == "__main__":
    calculate()
    a = input("press any keys to quit...")
