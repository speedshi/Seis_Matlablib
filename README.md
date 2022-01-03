# Seis_Matlablib
Matlab library which contains various functions and modules for seismology studies.

## Function list
### travel_times
**计算特定地震相位的旅行时**
1. tvtcalrt_ly: 计算层状介质中的直达波旅行时和take-off angle
2. tvtcalrt_homo: 计算均匀介质中的直达波旅行时和take-off angle

### model_build_show
**用于构建和显示模型的函数**
1. plotgeo: 用于显示三维层状介质模型，并在模型中同时显示震源、接收器阵列（包括表面阵列和垂直井接收阵列）的位置
2. plotmpsd: 用于绘制层状介质模型的速度、密度和衰减因子剖面
3. intpmodel3: 对输入的三维模型进行三维插值，得到符合期望(维度)大小的三维模型。注意输入三维模型的第一维为X，第二维为Y，第三维为Z
4. modelgeo_show: 根据输入参数，显示出相应模型和成像区域

### iofile
**用于输入输出文件的函数**
1. equihomo: 计算并输出层状介质的均方根速度
2. wrtasgeo: 输出震源，接收器位置信息
3. wrtmdf: 输出层状介质的速度、密度和衰减因子等信息，输出文件格式为fk软件输入格式
4. outputcatalogue: 输出地震目录，包括定位结果，地震时间和各台站到时
5. outputesgcsv: 按照ESG CSV 的格式输出地震信号
6. fdmodres: 读入fdmodeling合成的各道地震记录
7. rdfdmodres: 读入fdmodeling合成的各道地震记录，使用fdmodeling软件的模拟参数文件名和接收点文件名作为函数输入参数
8. wtraces2segy: 调用SegyMAT将合成地震记录输出为segy格式
9. read_seish5: 读入HDF5格式的地震数据，对HDF5文件的格式有固定的要求
10. read_stations: 读入IRIS text 格式的地震台站信息，获取台站的名称、位置信息
11. read_velocity: 读入均匀或层状介质的速度、层厚度信息
12. read_seissac: 读入SAC格式的地震数据
13. read_catalog: 读入IRIS text格式的地震目录
14. read_seis: main function to read seismic data of different formats (by calling different sub-functions)
15. read_staname: 读取HDF5数据文件中的台站名和台站数目
16. read_seismat: 读取MAT格式的地震数据
17. output_cataext1: 输出text格式catalog文件，适用于GMT画图


### seismic_modeling
**用于地震正演模拟的函数**
1. gsynwihl: 均匀介质或层状介质中的波形模拟。可以根据模型参数，自动选择调用均匀介质或层状介质的相关正演函数进行正演模拟。模型参数文件（震源文件、速度模型文件、接收器文件）的格式参照“file_format_description”。
2. rdmodelf: 读入均匀或层状介质模型参数文件(速度、密度、层厚度、衰减、起始深度)
3. rdreceiverf: 读入接收器参数文件
4. rdsourcef: 读入震源参数文件
5. plotmodel: 显示模型
6. example_build_input_files: 生成模拟所需的输入文件的示例代码
7. example_main_generate_waveforms: 调用模拟程序合成数据的示例代码
8. example_model.dat: 模型文件的示例
9. example_receiver.dat: 台站信息文件的示例
10. example_source.dat: 震源和模拟参数文件的示例
(1) homogeneous 均匀介质中的波场模拟
1. calninum: 计算Aki & Richards Eq 4.29 的解析解中近场的积分项，使用数值积分法
2. calnint:  计算Aki & Richards Equation 4.29 中的近场积分项，使用解析雷克子波
3. gsynwhomo: 利用Aki & Richards Eq 4.29 (解析解)计算均匀介质中的波场
4. gsynwhomo_rickerw: 利用Aki & Richards Eq 4.29 (解析解)计算均匀介质中的波场，使用解析雷克子波
5. homogreenf:  利用Aki & Richards Eq 4.29 (解析解)计算均匀介质中的波场,以脉冲为震源时间函数，计算出的是the Green's function
6. homogreenfne:  利用Aki & Richards Eq 4.29 (解析解)计算均匀介质中的近场波场,以脉冲为震源时间函数，计算出的是近场the Green's function (近场+中间场)
7. homogreenffa: 利用Aki & Richards Eq 4.29 (解析解)计算均匀介质中的远场波场,以脉冲为震源时间函数，计算出的是远场the Green's function
(2) layer 层状介质中的波场模拟
1. gsynwavefk/wavefk: 调用fk函数，实现层状介质中的反射率法正演模拟, 两个函数的输入参数略有不同。gsynwavefk函数更通用一些，更推荐使用该函数。



### general_math_func
**一般数学函数**
1. corrcoefn: 计算一个输入矩阵的n-维相关系数，n为矩阵的列数
2. corrcoefnv: 计算一个输入矩阵的多维相关系数
3. mycorrcoef: 计算输入矩阵的相关系数矩阵
4. mycovn: 计算输入矩阵的归一化协方差矩阵
5. mycroscorn: 计算输入矩阵的归一化互相关矩阵
6. my_kurtosis: 计算输入数据沿特定滑动时窗的kurtosis
7. my_stalta: 计算输入数据沿特定滑动时窗的STA/LTA
8. deltam: Dleta 函数
9. dnormlz: 将输入数据线性归一化到特定区间
10. gtwin: 生成特定窗函数的加权系数
11. intder: 计算输入数据的数值积分或微分
12. mtnorm: 将输入矩张量归一化
13. trdis2vel: 将地震数据由质点的位移分量转化为质点的速度分量
14. trvel2dis: 将地震数据由质点的速度分量转化为质点的位移分量
15. geod2cart: 将经纬度、海拔高度转化为笛卡尔坐标，使用wgs84Ellipsoid地理坐标系统
16. callyifdp: 计算层状介质的每一层界面（包括介质自由表面—起始深度）的深度
17. dnorm_mdn: 去除和获取数据的整体趋势，通过减去滑动时窗中的中位数来实现
18. datatransf: 对输入数据进行转换，如取绝对值，开方，取对数等


### display
**画图及显示相关函数**
1. disp_3dslice: 显示三个相互正交的剖面
2. migmaxplt: 显示一个输入4D 数据的最大值剖面和沿各维度的投影剖面，4D数据格式：T-X-Y-Z
3. disprs: 显示波剖面，即record section, 以震源接收器的水平距离为准排列道集
4. dispwfsc: 显示波剖面，类似record section，不同的是以震源接收器的直线距离为准排列道集，在记录的波形上标记P/S波到时
5. dispwfscn: 显示波剖面，类似record section，不同的是以震源接收器的直线距离为准排列道集，在每条记录的底线上标记P/S波到时
6. dispwfscn_2se: 同时显示两个波剖面，类似record section，不同的是以震源接收器的直线距离为准排列道集，在每条记录的底线上标记P/S波到时
7. particlemotion: 显示质点的振动轨迹
8. quiver3c: 显示三维矢量图（箭头），类似quiver3，矢量可以自由设置颜色
9. seisrsdisp: 显示地震剖面，按地震记录的顺序依次排列，每一道记录最大值归一化为1
10. seisrsdispk: 显示地震剖面，按地震记录的顺序依次排列，根据所有记录的最大值统一归一化为1
11. wigb: 显示地震波形记录
12. dispwflstk: 叠加并显示一定时窗内，波形的线性叠加结果
13. show_spectrogram: display the spectrogram of seismic data
14. ispectrogram/ispectrogram_1: display the seismogram and spectrogram of seismic data
15. dispwfscn_m:  显示地震波剖面,将台站名同时标注显示
16.  dispwfscn_mn: 显示地震波剖面,将台站名同时标注显示,台站按震源接收器的直线距离依次排列,间距为1
17. plot_evesta: 显示地震台站和地震事件的平面分布图

### seismic_location
**地震定位方法**  
(1) waveform_migration: 基于波形偏移的地震定位方法
1. stack_kernelf: 计算输入数据沿特定滑动时窗的特征函数
2. wavefmstk: 计算特征函数的叠加结果
3. mgrsprofdisp: 显示定位结果的XYZ剖面，并于地震目录中的结果对比
4. event_optm: 寻找偏移结果中的地震事件，采用时间、空间间隔的方式
5. extractevt: 寻找偏移结果中的地震事件，提取距离地震目录中的事件一定时间范围内的偏移最大值
6. locreson: 寻找偏移结果中的地震事件，提取在一定空间范围内持续一段时间的事件
7. findefmg: 寻找偏移结果中的地震事件，采用阈值和间隔时间的方式
8. gchkrs: 显示定位结果的记录剖面，帮助确认是否为明显的真实地震事件
9. profdisppw: 显示定位结果的XYZ剖面，并与地震目录中的地震事件一一对应
10. gpltlocrs: 在偏移记录上显示对应的定位结果（对应局部峰值）
11. mcm_genei: 读入各种数据，生成MCM Fortran程序所需的输入文件，并运行相应MCM程序
12. gene_soup: 生成偏移成像点的位置信息，并输出MCM需要的对应二进制文件
13. gene_traveltime: 生成旅行时表，并根据需要决定是否输出旅行时表二进制文件
14. gene_wavetime: 根据输入的波形数据生成MCM需要的波形二进制文件，并提取其对应的旅行时表并输出相应二进制文件
15. gene_migpara: 生成MCM所需的文本格式参数文件
16. runmcm_matlab_test: 根据输入的地震位置运行MCM matlab测试版本，会显示偏移剖面及记录剖面，用于判断偏移结果的好坏，可用于测试参数(如频率和时窗)的选择
17. waveform_migration_kernel: MCM偏移定位核心程序, use P+S phases
18. waveform_migration_kernel_x: MCM偏移定位核心程序, use only one phase
19. mcm_test_para.m: run MCM on a single position (soure location) to obtain the stacking trace to test the MCM parameters, such as frequency band, window size and seismic phases
20. get_earthquake: obtain the specified earthquake information from the catalog
21. mcm_test_freqband: test mcm results on different frequency bands
22. stkcorrcoef: calculate correlation coefficient matrix and stack the Ccs
23. waveforms_show: 显示定位结果的记录剖面，帮助确认是否为明显的真实地震事件,台站名和时间同时显示
24. gene_mcmifiles: 生成MCM Fortran 程序需要的各种参数文件（e.g. traveltime table, waveform file, soupos file, migpara）
25. detmst0: 生成搜索的orgin times 序列

### colormaps
**各种色标**
1. mycolor1.mat: 64*3, 蓝-黄-红
2. cmapmtrdneg.mat: 256*3, 蓝-黄
3. cmapmtdpos2: 256*3, 兰-黄-红
4. cmapmtrdpos: 64*3, 兰-黄-红
5. cmapmtrdp2: 256*3, 蓝-兰-黄-红
6. cmapmtrdp: 64*3, 蓝-兰-黄-红
7. cmapmtv: 64*3, 蓝-兰-黄-正红
8. cmapmtv2: 64*3, 蓝-兰-黄-正红, 兰黄占比增大
9. cmapmtv3: 64*3, 蓝-兰-黄-正红, 兰黄占比增大
10. cmapmtv4: 64*3, 蓝-兰-黄-正红, 兰黄占比增大

### downloads
**下载的各种函数和函数库**
1. borders: 显示世界各个国家的边界
2. color_map: 显示红蓝色标(地震剖面常用色标)
3. filter1: 对输入信号进行滤波
4. IPGP-sac_matlab-c67a67e: 对SAC文件进行读写
5. topotoolbox-master: 地形工具箱
6. XKCD_RGB: 获取不同颜色的RGB值
7. deg2utm: 将经纬度转化为UTM笛卡尔坐标
8. irisFetch: 链接IRIS，获取地震数据
9. freezeColors_v23_cbfreeze: 对不同子图使用不同的色标
10. segymat-1.6: 输入、输出和编辑segy格式的文件
11. fun_for_piero: clustering according to the input row and column indices of the upper triangular part of the correlation-coefficient matrix

### noise
**噪音有关函数**
1. addnoinsr: 按照噪信比（振幅比）向数据中加入指定噪音
2. pnoise: 向数据中加入一定信噪比的高斯随机噪音，信噪比以能量比表示
3. pnoisem: 向数据中加入一定噪信比的高斯随机噪音，噪信比以振幅比表示

### wavelet
**子波相关函数**
1. rickerw: 生成雷克子波，子波时延1.1/f+t0
2. rickerwd: 生成雷克子波导数（解析），子波时延1.1/f+t0
3. rickerwi: 生成雷克子波积分（数值），子波时延1.1/f+t0
4. waveldely: 将震源时间函数延迟
5. wavlintp: 将输入的震源时间函数插值加密，缩短时间间隔

### seismic source
**震源相关函数**
1. fgeom2mt: 由断层参数(走向、倾向、倾角)生成地震矩张量

### source radiation pattern
**震源辐射模式相关函数**
1. mtrdpfas: 画图相关函数，在画震源辐射模式图时，将坐标轴原点置于中心
2. mtrdpfax: 画图相关函数，在画震源辐射模式图时，将坐标轴原点置于中心；根据要画图像的数值自动选择合适的色标。
3. mtradiationv: 画远场P, S波的震源辐射模式，以矢量图(箭头)的形式展现
4.  mtradiationvbkv: 画远场P, S波的震源辐射模式，以矢量图(箭头)的形式展现，控制P波图中的矢量位置，使其图形更符合球形分布。
5. mtradiationb: 绘制远场P波的beach ball, 传统黑白beach ball, 三维球及使用stereographic projection的平面图（二维）
6. mtradiationb_prof: 绘制远场P波的beach ball, 蓝红色标, 三维球及其沿三个反坐标轴方向(-x, -y, -z)的三个二维平面视图
7. mtradiationp: 绘制远场P波，S波，SV波和SH波的震源辐射模式，三维图
8. mtradiationifps: 绘制中间场S波，远场P波，远场S波的震源辐射模式，三维图
9. mtradiationifps2: 绘制中间场S波，远场P波，远场S波在某个特定方位的震源辐射模式，相当与三维图的一个切片


### data_process
**一般性数据处理函数**
1. sltordotpsta: 按照与特定点的距离排列地震台站
2. stanamnum: 计算在输入的HFD5文件里，在特定station文件中的台站的数目和名称
3. wave_extract: extract waveforms along arrival times
4. ireplace_zeros: 使用随机的极小值（~eps）替换数据中的0值
5. decluster: 去除检测到的信号中连续存在的相同的事件，see code for detail
6. seisext: 根据输入的参数提取地震数据

### seismology
**地震学相关函数**
1. check_stations: 检测地震台站在一年中某个时间段是否有数据
2. catana_dist: 计算并显示地震目录中的地震时间距某点的距离
3. runFMF_mdata: 调用FMF的wrapper，注意其特定的数据输入模式
4. getpicks_fromNLLobs: 读取NonLinLoc格式的P/S 波的拾取文件
5. gene_FMFtemplate: 根据P/S波的拾取结果为FMF准备template
6. save_FMFdata: 保存FMF得到的检测结果，包括检测到的事件的事件和叠加相关系数

## Reference
If you use this package in your work, please cite the fellowing paper in your documents.  
```
Shi, P., Angus, D., Nowacki, A., Yuan, S., & Wang, Y. (2018). Microseismic full waveform modeling in anisotropic media with moment tensor implementation. Surveys in Geophysics, 39(4), 567-611.
```
```
@article{shi2018microseismic,
  title={Microseismic full waveform modeling in anisotropic media with moment tensor implementation},
  author={Shi, Peidong and Angus, Doug and Nowacki, Andy and Yuan, Sanyi and Wang, Yanyan},
  journal={Surveys in Geophysics},
  volume={39},
  number={4},
  pages={567--611},
  year={2018},
  publisher={Springer}
}
```

## Licence (GPLv3)
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. For more details, see in the license file.

## Contact information 
Copyright(C) 2021 Peidong Shi

Author: Peidong Shi

Email: peidong.shi@sed.ethz.ch or speedshi@hotmail.com
