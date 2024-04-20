from cmath import phase
import os
import sys
#*****************************************************************************
# Get the directory of the current script
main_directory = os.path.dirname(os.path.realpath(__file__))
parent_directory = os.path.dirname(main_directory)
parent_directory = os.path.dirname(parent_directory)

# Define the subdirectories you want to add to sys.path
subdirs = ['MaRGE', 'marcos_client']

# Add the subdirectories to sys.path
for subdir in subdirs:
    full_path = os.path.join(parent_directory, subdir)
    sys.path.append(full_path)

from typing import TypeAlias, Literal, Tuple, Dict, List
from enum import Enum
import numpy as np
import experiment as ex
import configs.hw_config as hw
import configs.units as units
from seq.mriBlankSeq import MRIBLANKSEQ as SeqBase

class 图像维度(Enum):
    读出 = 0
    相位 = 1
    选层 = 2

class B0Map(SeqBase):
    def __init__(self):
        super(B0Map, self).__init__()

        self.接收通道数 = 4
        self.激发脉冲时长 = 100  # μs
        self.激发脉冲幅度: float
        self.采样间隔: int

        self.nScans: int
        self.nPoints: Tuple[int, int, int]
        self.voxel: Tuple[float, float, float]
        self.axesOri: Tuple[int, int, int]
        self.larmorFreq: float
        self.flipAngle: float
        self.TR: float
        self.TE1: float
        self.TE2: float
        self.shimming: Tuple[float, float, float]
        self.gradRise: float
        self.gradRate: float
        
        参数类型: TypeAlias = str | int | float | List[int | float]
        标签类型: TypeAlias = Literal['IM'] | Literal['RF'] | Literal['SEQ'] | Literal['OTH']
        序列参数: Dict[str, Tuple[str, 参数类型, 标签类型, float, str]] = {
            'seqName': ['B0MapInfo', 'B0Map'],
            'nScans': ['平均次数', 1, 'IM'],
            'nPoints': ['采样点数(读出,相位,选层)', [64, 64, 1], 'IM'],
            'voxel': ['体素 (mm)', [1.0, 1.0, 1.0], 'IM', units.mm],
            'axesOri': ['Axes orientation', [0, 1, 2], 'IM'],
            'larmorFreq': ['Larmor frequency (MHz)', hw.larmorFreq, 'RF'],
            'flipAngle': ['翻转角 (°)', 90, 'RF'],
            'TR': ['Repertition time (ms)', 1000, 'SEQ', units.ms],
            'TE1': ['Echo time 1 (ms)', 1.0, 'SEQ', units.ms],
            'TE2': ['Echo time 2 (ms)', 2.0, 'SEQ', units.ms],
            'shimming': ['Shimming (*1e4)', [0.0, 0.0, 0.0], 'OTH', units.sh],
            'gradRise': ['Gradient rising (μs)', 50, 'OTH', units.us],
            'gradRate': ['Gradient sampling rate (MHz)', 0.2, 'OTH'],
        }
        for key, value in 序列参数.items():
            self.addParameter(
                key, 
                string=value[0], 
                val=value[1],
                field=value[2] if len(value) > 2 else None,
                units=value[3] if len(value) > 3 else True
            )


    def sequenceInfo(self):
        print("""
            一次激发采集两个梯度回波用于计算B0场。
        """.lstrip())


    def sequenceTime(self) -> float:
        self.sequenceAtributes()
        return self.nScans * self.相位编码步数 * self.层数 * self.TR / 60


    def sequenceRun(self, plotSeq=False, demo=False):
        self.统一时间单位()
        if self.读出时长 <= 0 or self.回波间隔 < 0:
            print('TE is too short')
            return False
        带宽 = self.读出点数 / self.读出时长
        print('带宽: %.2f kHz' % (带宽 * 1e3))

        try:
            self.采样间隔 = self.dummyPulse(1 / 带宽)
            self.原始数据 = []

            for i in range(self.nScans):
                self.expt.__del__()
                self.expt = ex.Experiment(
                    self.larmorFreq, 
                    self.采样间隔, 
                    grad_max_update_rate=self.gradRate
                )
                self.写序列()
                if self.floDict2Exp():
                    print('正在扫描第%d次' % (i + 1))
                    rxd, _ = self.expt.run()
                    self.原始数据 += [list(rxd.values())]#np.array([rxd['rx%d' % j] for j in range(self.接收通道数)])
                else:
                    print('sequence compile error')
                    return False
        except Exception as e:
            print(e)
            return False
        finally:
            self.expt.__del__()
        return True
        
    
    def sequenceAnalysis(self, mode: Literal['Standalone']=None):
        def kSpace(raw: np.ndarray) -> np.ndarray:
            dualEcho = raw.reshape(-1, self.读出点数)
            for i in range(dualEcho.shape[0]):
                dualEcho[i] *= np.exp(-1j * self.共振相位(i))
            dualEcho = dualEcho.reshape(self.层数, self.相位编码步数, -1)
            echo1 = dualEcho[:, :, :self.读出点数]
            echo2 = dualEcho[:, :, self.读出点数:]
            return np.array([echo1, echo2])

        self.平均_通道_回波_选层_相位_读出 = [[
            kSpace(rxd) for rxd in self.原始数据[i]
        ] for i in range(self.nScans)]
        ksp5d = np.mean(self.平均_通道_回波_选层_相位_读出, axis=0)
        img5d = self.mapVals['img5d'] = [[
            self.重建(ksp5d[i, j]) for j in range(2)
        ] for i in range(self.接收通道数)]

        if self.axesOri[图像维度.选层.value] < 0 and self.axesOri[图像维度.相位.value] < 0:
            self.output = [{
                'widget': 'curve',
                'xLabel': 'time (us)',
                'xData': np.linspace(0, self.采样间隔 * self.读出点数, self.读出点数),
                'yLabel': 'signal',
                'yData': [curve.real, curve.imag, np.abs(curve), self.相位解缠绕(curve)],
                'title': 'in-phase' if i == 0 else 'out-of-phase',
                'legend': ['real', 'imaginary', 'abs', 'phase'],
                'row': 0,
                'col': i
            } for i, curve in enumerate(ksp5d[0, :, 0, 0])]
        else:
            from skimage.restoration import unwrap_phase
            from skimage.transform import hough_ellipse
            from skimage.feature import canny
            from skimage.draw import ellipse, ellipse_perimeter
            
            img3d1 = img5d[0][0]
            img3d2 = img5d[0][1][..., ::-1]  # reverse readout for TE2

            abs3d1 = np.abs(img3d1)
            abs3d2 = np.abs(img3d2)
            bin3d = [canny(img, 2, 0.001, 0.3, use_quantiles=True) for img in abs3d1]  # 第一个回波比较清晰
            phase3d1 = []
            phase3d2 = []
            
            for i in range(self.层数):
                ellipses = hough_ellipse(bin3d[i], min_size=12, threshold=12)
                try:
                    过滤 = [e['a'] > 6 and e['b'] > 6 for e in ellipses]
                    ellipses = ellipses[过滤]
                    ellipses.sort(order='accumulator')
                    best = list(ellipses[-1])
                    yc, xc, a, b = (int(round(i)) for i in best[1:5])
                    rot = best[5]
                    rr, cc = ellipse_perimeter(yc, xc, a, b, rot, mask.shape)
                    abs3d2[i][rr, cc] = abs3d2[i].max() * 2

                    mask = np.zeros_like(bin3d[i])
                    mask[ellipse(yc, xc, a, b, bin3d[i].shape, rot)] = 1
                except Exception as e:
                    print('no ellipse found in layer %d' % i)
                    print(e)
                    mask = np.ones_like(bin3d[i])
                    continue
                finally:
                    phase3d1 += [unwrap_phase(np.angle(img3d1[i]) * mask)]
                    phase3d2 += [unwrap_phase(np.angle(img3d2[i]) * mask)]
                    
            Δɸ = np.array(phase3d2) - np.array(phase3d1)
            b0 = Δɸ / self.ΔTE
            self.output = [{
                'widget': 'image',
                'data': np.concatenate((b0, bin3d * b0.max()), axis=1),
                'title': 'B0 map (Hz)',
                'xLabel': 'x',
                'yLabel': 'y',
                'row': 0,
                'col': 0
            }, {
                'widget': 'image',
                'data': np.concatenate((abs3d1, abs3d2), axis=1),
                'title': 'image',
                'xLabel': 'x',
                'yLabel': 'y',
                'row': 0,
                'col': 1
            }, {
                'widget': 'image',
                'data': np.concatenate((phase3d1, phase3d2), axis=1),
                'title': 'phase',
                'xLabel': 'kx',
                'yLabel': 'ky',
                'row': 1,
                'col': 0
            }, {
                'widget': 'image',
                'data': np.concatenate(np.abs(ksp5d[0][:]), axis=1),
                'title': 'ksp',
                'xLabel': 'kx',
                'yLabel': 'ky',
                'row': 1,
                'col': 1
            }]
        self.保存()
        return self.output
    
###############################################################################
    
    @property
    def 读出点数(self) -> int:
        return self.nPoints[图像维度.读出.value]

    @property
    def 相位编码步数(self) -> int:
        return self.nPoints[图像维度.相位.value]

    @property
    def 层数(self) -> int:
        return self.nPoints[图像维度.选层.value]

    @property
    def 激发脉冲幅度(self) -> float:
        return self.翻转角换算幅度(self.激发脉冲时长)
    
    @property
    def 读出时长(self) -> float:
        return self.TE1 - self.激发脉冲时长/2 - 2*self.gradRise
    
    @property
    def 编码时长(self) -> float:
        return self.读出时长/2 - self.gradRise
    
    @property
    def ΔTE(self) -> float:
        return self.TE2 - self.TE1

    @property
    def 回波间隔(self) -> float:
        return self.ΔTE - 2*self.gradRise - self.读出时长
    
    
    def 相位解缠绕(self, curve: np.ndarray) -> np.ndarray:
        unwrapped = np.unwrap(np.angle(curve))
        return unwrapped / np.max(np.abs(unwrapped))


    def 重建(self, data: np.ndarray):
        return np.fft.ifftshift(np.fft.ifftn(np.fft.ifftshift(data)))


    def 翻转角换算幅度(self, 脉冲时长=100) -> float:
        return self.flipAngle / 180 / np.pi / 脉冲时长 / hw.b1Efficiency


    def 梯度编码幅值(self, axes: 图像维度) -> float:
        return hw.gFactor[self.axesOri[axes.value]] / hw.gammaB / self.voxel[axes.value] / (self.编码时长 + self.gradRise) * 1e6


    def dummyPulse(self, 采样间隔: float) -> float:
        """
        初始化梯度板，平衡磁化强度，返回采样间隔（有副作用）
        """
        self.iniSequence(20, self.shimming)
        self.rfRecPulse(1e4, self.激发脉冲时长, self.激发脉冲幅度)
        self.rxGate(1e4 + self.激发脉冲时长, self.读出时长)
        self.endSequence(self.TR + 2e4)
        self.expt = ex.Experiment(
            self.larmorFreq, 
            采样间隔, 
            grad_max_update_rate=self.gradRate,
            init_gpa=True
        )
        if self.floDict2Exp():
            _, msg = self.expt.run()
            print(msg)
        else:
            print('dummy pulse failed')
        return self.expt.get_rx_ts()[0] #* hw.oversamplingFactor


    def 统一时间单位(self):
        self.sequenceAtributes()
        self.gradRise /= units.us
        self.TR /= units.us
        self.TE1 /= units.us
        self.TE2 /= units.us


    def 保存(self):
        self.mapVals['raw'] = self.原始数据
        self.mapVals['ksp'] = self.平均_通道_回波_选层_相位_读出
        self.mapVals['out'] = self.output
        self.mapVals['bw'] = self.读出点数 / self.读出时长 * 1e6
        self.saveRawData()


    def 共振相位(self, idx: int) -> float:
        # return np.pi/2 * (1 + (-1)**idx)
        return 117 * np.pi / 180 * idx


    def 写序列(self):
        
        def 射频脉冲(t0: float, idx: int):
            self.rfRecPulse(t0, self.激发脉冲时长, self.激发脉冲幅度, self.共振相位(idx))
            return hw.blkTime + self.激发脉冲时长
        
        def 接收门控(t0: float):
            for ch in range(self.接收通道数):
                self.rxGate(t0, self.采样间隔 * self.读出点数, ch)
            return self.读出时长
        
        def 梯形梯度(t0: float, 平台时长: float, 幅度: float, 通道: int):
            if 通道 >= 0:
                self.gradTrap(
                    t0, 
                    self.gradRise, 
                    平台时长, 
                    幅度, 
                    hw.grad_steps,
                    通道,
                    self.shimming
                )
            return self.gradRise*2 + 平台时长
        
        kxAmp = 1e6 * hw.gFactor[self.axesOri[图像维度.读出.value]] / hw.gammaB / self.voxel[图像维度.读出.value] / self.读出时长
        kyAmp = self.梯度编码幅值(图像维度.相位)
        相位编码 = np.linspace(-kyAmp, kyAmp, self.相位编码步数)
        kzAmp = self.梯度编码幅值(图像维度.选层)
        选层编码 = np.linspace(-kzAmp, kzAmp, self.层数)
        print('gx: %f, gy: %f, gz: %f' % (kxAmp, kyAmp, kzAmp))
        
        self.iniSequence(20, self.shimming)
        t0 = 1e4

        for i in range(self.层数):
            for j in range(self.相位编码步数):
                t0 += self.TR
                t = t0 + 射频脉冲(t0, i*self.相位编码步数 + j)

                # 相位编码
                梯形梯度(
                    t,
                    self.编码时长,
                    选层编码[i],
                    self.axesOri[图像维度.选层.value]
                )
                梯形梯度(
                    t + 1,
                    self.编码时长,
                    相位编码[j],
                    self.axesOri[图像维度.相位.value]
                )
                
                # 频率编码
                t += 梯形梯度(
                    t, 
                    self.编码时长, 
                    kxAmp, 
                    self.axesOri[图像维度.读出.value]
                )

                # 同相回波 (in-phase)
                接收门控(t + self.gradRise)
                t += 梯形梯度(
                    t, 
                    self.读出时长, 
                    -kxAmp, 
                    self.axesOri[图像维度.读出.value]
                )

                # 异相回波 (out-of-phase)
                接收门控(t + self.gradRise + self.回波间隔)
                t += 梯形梯度(
                    t + self.回波间隔, 
                    self.读出时长, 
                    kxAmp, 
                    self.axesOri[图像维度.读出.value]
                )

        self.endSequence(t0 + self.TR + 6e4)
