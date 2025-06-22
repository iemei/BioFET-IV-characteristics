import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy as np
from scipy.stats import linregress

data = pd.read_csv(r"C:\Users\iemei\PycharmProjects\PythonProject\Ids-Vbg_Vjg=-1.5Vto0.5V.csv")
print(data.head())
voltage = data[" Vbg"]
current1 = abs(data[" Id_Vjg=-1.5V"])
current2 = abs(data[" Id_Vjg=-1.25V"])
current3 = abs(data[" Id_Vjg=-1V"])
current4 = abs(data[" Id_Vjg=-0.75V"])
current5 = abs(data[" Id_Vjg=-0.5V"])
current6 = abs(data[" Id_Vjg=-0.25V"])
current7 = abs(data[" Id_Vjg=0V"])

current1_smooth = savgol_filter(current1, window_length=11, polyorder=1)
current2_smooth = savgol_filter(current2, window_length=45, polyorder=1)
current3_smooth = savgol_filter(current3, window_length=25, polyorder=1)
current4_smooth = savgol_filter(current4, window_length=25, polyorder=1)
current5_smooth = savgol_filter(current5, window_length=25, polyorder=1)
current6_smooth = savgol_filter(current6, window_length=25, polyorder=1)
current7_smooth = savgol_filter(current7, window_length=25, polyorder=1)

plt.plot(voltage, current1, alpha=0.5, color='orange')
plt.plot(voltage, current1_smooth, label='Smoothed Vjg = -1.5V', linewidth=2, color='orange', linestyle='--')
plt.plot(voltage, current2, alpha=0.5, color='blue')
plt.plot(voltage, current2_smooth, label='Smoothed Vjg = -1.25V', linewidth=2, color='blue', linestyle='--')
plt.plot(voltage, current3, alpha=0.5, color='green')
plt.plot(voltage, current3_smooth, label='Smoothed Vjg = -1.0V', linewidth=2, color='green', linestyle='--')
plt.plot(voltage, current4, alpha=0.5, color='red')
plt.plot(voltage, current4_smooth, label='Smoothed Vjg = -0.75V', linewidth=2, color='red', linestyle='--')
plt.plot(voltage, current5, alpha=0.5, color='black')
plt.plot(voltage, current5_smooth, label='Smoothed Vjg = -0.5V', linewidth=2, color='black', linestyle='--')
plt.plot(voltage, current6, alpha=0.5, color='purple')
plt.plot(voltage, current6_smooth, label='Smoothed Vjg = -0.25V', linewidth=2, color='purple', linestyle='--')
plt.plot(voltage, current7, alpha=0.5, color='magenta')
plt.plot(voltage, current7_smooth, label='Smoothed Vjg = 0V', linewidth=2, color='magenta', linestyle='--')

plt.xlabel("Backgate Voltage Vbg (V)")
plt.ylabel("Drain Current Ids (A)")
plt.yscale("log")
plt.legend()
plt.grid(True)
plt.title("Smoothed I–V Curve")
plt.show()

#Threshold voltage extraction

sqrt_ids1 = np.sqrt(current1_smooth)
sqrt_ids2 = np.sqrt(current2_smooth)
sqrt_ids3 = np.sqrt(current3_smooth)
sqrt_ids4 = np.sqrt(current4_smooth)
sqrt_ids5 = np.sqrt(current5_smooth)
sqrt_ids6 = np.sqrt(current6_smooth)
sqrt_ids7 = np.sqrt(current7_smooth)

fit_region = (voltage > -9) & (voltage < -8)

slope1, intercept1, r_value1, p_value1, std_err1 = linregress(voltage[fit_region], sqrt_ids1[fit_region])
vth1 = -intercept1 / slope1
print(f"Extracted Vth = {vth1:0.3f} V")

slope2, intercept2, r_value2, p_value2, std_err2 = linregress(voltage[fit_region], sqrt_ids2[fit_region])
vth2 = -intercept2 / slope2
print(f"Extracted Vth = {vth2:0.3f} V")

slope3, intercept3, r_value3, p_value3, std_err3 = linregress(voltage[fit_region], sqrt_ids3[fit_region])
vth3 = -intercept3 / slope3
print(f"Extracted Vth = {vth3:0.3f} V")

slope4, intercept4, r_value4, p_value4, std_err4 = linregress(voltage[fit_region], sqrt_ids4[fit_region])
vth4 = -intercept4 / slope4
print(f"Extracted Vth = {vth4:0.3f} V")

slope5, intercept5, r_value5, p_value5, std_err5 = linregress(voltage[fit_region], sqrt_ids5[fit_region])
vth5 = -intercept5 / slope5
print(f"Extracted Vth = {vth5:0.3f} V")

slope6, intercept6, r_value6, p_value6, std_err6 = linregress(voltage[fit_region], sqrt_ids6[fit_region])
vth6 = -intercept6 / slope6
print(f"Extracted Vth = {vth6:0.3f} V")

slope7, intercept7, r_value7, p_value7, std_err7 = linregress(voltage[fit_region], sqrt_ids7[fit_region])
vth7 = -intercept7 / slope7
print(f"Extracted Vth = {vth7:0.3f} V")

plt.plot(voltage, sqrt_ids1, label="Vjg=-1.5V", color='orange')
plt.plot(voltage[fit_region], slope1 * voltage[fit_region] + intercept1, 'r--')
plt.plot(voltage, sqrt_ids2, label="Vjg=-1.25V", color='blue')
plt.plot(voltage[fit_region], slope2 * voltage[fit_region] + intercept2, 'r--')
plt.plot(voltage, sqrt_ids3, label="Vjg=-1.0V", color='green')
plt.plot(voltage[fit_region], slope3 * voltage[fit_region] + intercept3, 'r--')
plt.plot(voltage, sqrt_ids4, label="Vjg=-0.75V", color='red')
plt.plot(voltage[fit_region], slope4 * voltage[fit_region] + intercept4, 'r--')
plt.plot(voltage, sqrt_ids5, label="Vjg=-0.5V", color='black')
plt.plot(voltage[fit_region], slope5 * voltage[fit_region] + intercept5, 'r--')
plt.plot(voltage, sqrt_ids6, label="Vjg=-0.25V", color='purple')
plt.plot(voltage[fit_region], slope6 * voltage[fit_region] + intercept6, 'r--')
plt.plot(voltage, sqrt_ids7, label="Vjg=0V", color='magenta')
plt.plot(voltage[fit_region], slope7 * voltage[fit_region] + intercept7, 'r--')

plt.yscale("log")
plt.axvline(vth1, color='orange', linestyle=':', label=f"Vth ≈ {vth1:0.3f} V")
plt.axvline(vth2, color='blue', linestyle=':', label=f"Vth ≈ {vth2:0.3f} V")
plt.axvline(vth3, color='green', linestyle=':', label=f"Vth ≈ {vth3:0.3f} V")
plt.axvline(vth4, color='red', linestyle=':', label=f"Vth ≈ {vth4:0.3f} V")
plt.axvline(vth5, color='black', linestyle=':', label=f"Vth ≈ {vth5:0.3f} V")
plt.axvline(vth6, color='purple', linestyle=':', label=f"Vth ≈ {vth6:0.3f} V")
plt.axvline(vth7, color='magenta', linestyle=':', label=f"Vth ≈ {vth7:0.3f} V")
plt.xlabel("Backgate Voltage Vbg (V)")
plt.ylabel("√Ids")
plt.legend()
plt.grid(True)
plt.title("Threshold Voltage Extraction (√Id method)")
plt.show()

gm1 = np.gradient(current1_smooth, voltage)
gm1_smooth = savgol_filter(gm1, window_length=75, polyorder=1)
gm2 = np.gradient(current2_smooth, voltage)
gm2_smooth = savgol_filter(gm2, window_length=49, polyorder=1)
gm3 = np.gradient(current3_smooth, voltage)
gm3_smooth = savgol_filter(gm3, window_length=57, polyorder=1)
gm4 = np.gradient(current4_smooth, voltage)
gm4_smooth = savgol_filter(gm4, window_length=53, polyorder=1)
gm5 = np.gradient(current5_smooth, voltage)
gm5_smooth = savgol_filter(gm5, window_length=51, polyorder=1)
gm6 = np.gradient(current6_smooth, voltage)
gm6_smooth = savgol_filter(gm6, window_length=41, polyorder=1)
gm7 = np.gradient(current7_smooth, voltage)
gm7_smooth = savgol_filter(gm7, window_length=35, polyorder=1)

gm_max_index1 = np.argmax(gm1_smooth)
vth_gm1 = voltage.iloc[gm_max_index1]
gm_max_index2 = np.argmax(gm2_smooth)
vth_gm2 = voltage.iloc[gm_max_index2]
gm_max_index3 = np.argmax(gm3_smooth)
vth_gm3 = voltage.iloc[gm_max_index3]
gm_max_index4 = np.argmax(gm4_smooth)
vth_gm4 = voltage.iloc[gm_max_index4]
gm_max_index5 = np.argmax(gm5_smooth)
vth_gm5 = voltage.iloc[gm_max_index5]
gm_max_index6 = np.argmax(gm6_smooth)
vth_gm6 = voltage.iloc[gm_max_index6]
gm_max_index7 = np.argmax(gm7_smooth)
vth_gm7 = voltage.iloc[gm_max_index7]

plt.figure(figsize=(8, 5))
plt.plot(voltage, gm1_smooth, label='Vjg=-1.5V', color='orange')
plt.plot(voltage, gm2_smooth, label='Vjg=-1.25V', color='blue')
plt.plot(voltage, gm3_smooth, label='Vjg=-1.0V', color='green')
plt.plot(voltage, gm4_smooth, label='Vjg=-0.75V', color='red')
plt.plot(voltage, gm5_smooth, label='Vjg=-0.5V', color='black')
plt.plot(voltage, gm6_smooth, label='Vjg=-0.25V', color='purple')
plt.plot(voltage, gm7_smooth, label='Vjg=0V', color='magenta')
plt.axvline(vth_gm1, color='orange', linestyle='--', label=f'V_th ≈ {vth_gm1:.3f} V')
plt.axvline(vth_gm2, color='blue', linestyle='--', label=f'V_th ≈ {vth_gm2:.3f} V')
plt.axvline(vth_gm3, color='green', linestyle='--', label=f'V_th ≈ {vth_gm3:.3f} V')
plt.axvline(vth_gm4, color='red', linestyle='--', label=f'V_th ≈ {vth_gm4:.3f} V')
plt.axvline(vth_gm5, color='black', linestyle='--', label=f'V_th ≈ {vth_gm5:.3f} V')
plt.axvline(vth_gm6, color='purple', linestyle='--', label=f'V_th ≈ {vth_gm6:.3f} V')
plt.axvline(vth_gm7, color='magenta', linestyle='--', label=f'V_th ≈ {vth_gm7:.3f} V')

plt.xlabel("Backgate Voltage Vbg (V)")
plt.ylabel("gₘ = dI_D/dV_G (A/V)")
plt.title("Threshold Voltage via Maximum Transconductance Method")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()
print(f"Extracted Threshold Voltage (Max gm method): {vth_gm1:0.3f} V")
print(f"Extracted Threshold Voltage (Max gm method): {vth_gm2:0.3f} V")
print(f"Extracted Threshold Voltage (Max gm method): {vth_gm3:0.3f} V")
print(f"Extracted Threshold Voltage (Max gm method): {vth_gm4:0.3f} V")
print(f"Extracted Threshold Voltage (Max gm method): {vth_gm5:0.3f} V")
print(f"Extracted Threshold Voltage (Max gm method): {vth_gm6:0.3f} V")
print(f"Extracted Threshold Voltage (Max gm method): {vth_gm7:0.3f} V")

dgm1 = np.gradient(gm1_smooth, voltage)
dgm1_smooth = savgol_filter(dgm1, window_length=95, polyorder=1)
dgm1_max_index = np.argmax(dgm1_smooth)
#zero_crossings = np.where(np.diff(np.sign(dgm1_smooth)))[0]
vth_second_deriv1 = voltage.iloc[dgm1_max_index]

dgm2 = np.gradient(gm2_smooth, voltage)
dgm2_smooth = savgol_filter(dgm2, window_length=91, polyorder=1)
dgm2_max_index = np.argmax(dgm2_smooth)
#zero_crossings = np.where(np.diff(np.sign(dgm2_smooth)))[0]
vth_second_deriv2 = voltage.iloc[dgm2_max_index]

dgm3 = np.gradient(gm3_smooth, voltage)
dgm3_smooth = savgol_filter(dgm3, window_length=87, polyorder=1)
dgm3_max_index = np.argmax(dgm3_smooth)
#zero_crossings = np.where(np.diff(np.sign(dgm3_smooth)))[0]
vth_second_deriv3 = voltage.iloc[dgm3_max_index]

dgm4 = np.gradient(gm4_smooth, voltage)
dgm4_smooth = savgol_filter(dgm4, window_length=79, polyorder=1)
dgm4_max_index = np.argmax(dgm4_smooth)
#zero_crossings = np.where(np.diff(np.sign(dgm4_smooth)))[0]
vth_second_deriv4 = voltage.iloc[dgm4_max_index]

dgm5 = np.gradient(gm5_smooth, voltage)
dgm5_smooth = savgol_filter(dgm5, window_length=75, polyorder=1)
dgm5_max_index = np.argmax(dgm5_smooth)
#zero_crossings = np.where(np.diff(np.sign(dgm5_smooth)))[0]
vth_second_deriv5 = voltage.iloc[dgm5_max_index]

dgm6 = np.gradient(gm6_smooth, voltage)
dgm6_smooth = savgol_filter(dgm6, window_length=75, polyorder=1)
dgm6_max_index = np.argmax(dgm6_smooth)
#zero_crossings = np.where(np.diff(np.sign(dgm6_smooth)))[0]
vth_second_deriv6 = voltage.iloc[dgm6_max_index]

dgm7 = np.gradient(gm7_smooth, voltage)
dgm7_smooth = savgol_filter(dgm7, window_length=75, polyorder=1)
dgm7_max_index = np.argmax(dgm7_smooth)
#zero_crossings = np.where(np.diff(np.sign(dgm7_smooth)))[0]
vth_second_deriv7 = voltage.iloc[dgm7_max_index]


plt.figure(figsize=(10, 6))
plt.plot(voltage, dgm1_smooth, label="Vjg = -1.5V", color='orange')
plt.plot(voltage, dgm2_smooth, label="Vjg = -1.25V", color='blue')
plt.plot(voltage, dgm3_smooth, label="Vjg = -1.0V", color='green')
plt.plot(voltage, dgm4_smooth, label="Vjg = -0.75V", color='red')
plt.plot(voltage, dgm5_smooth, label="Vjg = -0.5V", color='black')
plt.plot(voltage, dgm6_smooth, label="Vjg = -0.25V", color='purple')
plt.plot(voltage, dgm7_smooth, label="Vjg = 0V", color='magenta')

if vth_second_deriv1:
    plt.axvline(vth_second_deriv1, color='orange', linestyle='--', label=f"V_th ≈ {vth_second_deriv1:.3f} V")
if vth_second_deriv2:
    plt.axvline(vth_second_deriv2, color='blue', linestyle='--', label=f"V_th ≈ {vth_second_deriv2:.3f} V")
if vth_second_deriv3:
    plt.axvline(vth_second_deriv3, color='green', linestyle='--', label=f"V_th ≈ {vth_second_deriv3:.3f} V")
if vth_second_deriv4:
    plt.axvline(vth_second_deriv4, color='red', linestyle='--', label=f"V_th ≈ {vth_second_deriv4:.3f} V")
if vth_second_deriv5:
    plt.axvline(vth_second_deriv5, color='black', linestyle='--', label=f"V_th ≈ {vth_second_deriv5:.3f} V")
if vth_second_deriv6:
    plt.axvline(vth_second_deriv6, color='purple', linestyle='--', label=f"V_th ≈ {vth_second_deriv6:.3f} V")
if vth_second_deriv7:
    plt.axvline(vth_second_deriv7, color='magenta', linestyle='--', label=f"V_th ≈ {vth_second_deriv7:.3f} V")

plt.xlabel("Backgate Voltage Vbg (V)")
plt.ylabel("d²Ids / dVbg²")
plt.title("Threshold Voltage via Second Derivative Method")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

if vth_second_deriv1:
    print(f"Extracted Threshold Voltage (2nd derivative method): {vth_second_deriv1:.3f} V")
else:
    print("Could not detect zero crossing in second derivative.")
if vth_second_deriv2:
    print(f"Extracted Threshold Voltage (2nd derivative method): {vth_second_deriv2:.3f} V")
else:
    print("Could not detect zero crossing in second derivative.")
if vth_second_deriv3:
    print(f"Extracted Threshold Voltage (2nd derivative method): {vth_second_deriv3:.3f} V")
else:
    print("Could not detect zero crossing in second derivative.")
if vth_second_deriv4:
    print(f"Extracted Threshold Voltage (2nd derivative method): {vth_second_deriv4:.3f} V")
else:
    print("Could not detect zero crossing in second derivative.")
if vth_second_deriv5:
    print(f"Extracted Threshold Voltage (2nd derivative method): {vth_second_deriv5:.3f} V")
else:
    print("Could not detect zero crossing in second derivative.")
if vth_second_deriv6:
    print(f"Extracted Threshold Voltage (2nd derivative method): {vth_second_deriv6:.3f} V")
else:
    print("Could not detect zero crossing in second derivative.")
if vth_second_deriv7:
    print(f"Extracted Threshold Voltage (2nd derivative method): {vth_second_deriv7:.3f} V")
else:
    print("Could not detect zero crossing in second derivative.")

# SS extraction

log_current1 = np.log10(current1_smooth)
log_current2 = np.log10(current2_smooth)
log_current3 = np.log10(current3_smooth)
log_current4 = np.log10(current4_smooth)
log_current5 = np.log10(current5_smooth)
log_current6 = np.log10(current6_smooth)
log_current7 = np.log10(current7_smooth)

subregion = (voltage > -9) & (voltage < -8)
slope11, intercept11, r_value11, _, _ = linregress(voltage[subregion], log_current1[subregion])
slope12, intercept12, r_value12, _, _ = linregress(voltage[subregion], log_current2[subregion])
slope13, intercept13, r_value13, _, _ = linregress(voltage[subregion], log_current3[subregion])
slope14, intercept14, r_value14, _, _ = linregress(voltage[subregion], log_current4[subregion])
slope15, intercept15, r_value15, _, _ = linregress(voltage[subregion], log_current5[subregion])
slope16, intercept16, r_value16, _, _ = linregress(voltage[subregion], log_current6[subregion])
slope17, intercept17, r_value17, _, _ = linregress(voltage[subregion], log_current7[subregion])
SS1 = 1 / slope11 * 1000 * np.log(10)
r_squared1 = r_value11**2
print(f"Subthreshold Swing (SS) ≈ {SS1:.2f} mV/dec")
print(f"Subthreshold Linearity (R^2) ≈ {r_squared1:.4f}")
SS2 = 1 / slope12 * 1000 * np.log(10)
r_squared2 = r_value12**2
print(f"Subthreshold Swing (SS) ≈ {SS2:.2f} mV/dec")
print(f"Subthreshold Linearity (R^2) ≈ {r_squared2:.4f}")
SS3 = 1 / slope13 * 1000 * np.log(10)
r_squared3 = r_value13**2
print(f"Subthreshold Swing (SS) ≈ {SS3:.2f} mV/dec")
print(f"Subthreshold Linearity (R^2) ≈ {r_squared3:.4f}")
SS4 = 1 / slope14 * 1000 * np.log(10)
r_squared4 = r_value14**2
print(f"Subthreshold Swing (SS) ≈ {SS4:.2f} mV/dec")
print(f"Subthreshold Linearity (R^2) ≈ {r_squared4:.4f}")
SS5 = 1 / slope15 * 1000 * np.log(10)
r_squared5 = r_value15**2
print(f"Subthreshold Swing (SS) ≈ {SS5:.2f} mV/dec")
print(f"Subthreshold Linearity (R^2) ≈ {r_squared5:.4f}")
SS6 = 1 / slope16 * 1000 * np.log(10)
r_squared6 = r_value16**2
print(f"Subthreshold Swing (SS) ≈ {SS6:.2f} mV/dec")
print(f"Subthreshold Linearity (R^2) ≈ {r_squared6:.4f}")
SS7 = 1 / slope17 * 1000 * np.log(10)
r_squared7 = r_value17**2
print(f"Subthreshold Swing (SS) ≈ {SS7:.2f} mV/dec")
print(f"Subthreshold Linearity (R^2) ≈ {r_squared7:.4f}")

plt.plot(voltage, log_current1, label='log10(Ids)', color='orange')
plt.plot(voltage[subregion], slope11 * voltage[subregion] + intercept11, 'r--', label='Subthreshold Fit1')
plt.plot(voltage, log_current2, label='log10(Ids)', color='blue')
plt.plot(voltage[subregion], slope12 * voltage[subregion] + intercept12, 'r--', label='Subthreshold Fit2')
plt.plot(voltage, log_current3, label='log10(Ids)',color='green')
plt.plot(voltage[subregion], slope13 * voltage[subregion] + intercept13, 'r--', label='Subthreshold Fit3')
plt.plot(voltage, log_current4, label='log10(Ids)', color='red')
plt.plot(voltage[subregion], slope14 * voltage[subregion] + intercept14, 'r--', label='Subthreshold Fit4')
plt.plot(voltage, log_current5, label='log10(Ids)', color='black')
plt.plot(voltage[subregion], slope15 * voltage[subregion] + intercept15, 'r--', label='Subthreshold Fit5')
plt.plot(voltage, log_current6, label='log10(Ids)', color='purple')
plt.plot(voltage[subregion], slope16 * voltage[subregion] + intercept16, 'r--', label='Subthreshold Fit6')
plt.plot(voltage, log_current7, label='log10(Ids)', color='magenta')
plt.plot(voltage[subregion], slope17 * voltage[subregion] + intercept17, 'r--', label='Subthreshold Fit7')


plt.xlabel("Vgs (V)")
plt.ylabel("log10(Ids)")
plt.title("Subthreshold Region Fit")
plt.grid(True)
plt.legend()
plt.show()

# on-off ratio
I_on1 = current1.max()
I_off1 = current1.min()
on_off_ratio1 = I_on1/I_off1
print(f"On Current I_on ≈ {I_on1:.3e} A ")
print(f"On Current I_off ≈ {I_off1:.3e} A")
print(f"Current on-off ratio ≈ {on_off_ratio1:.2e} ")
I_on2 = current2.max()
I_off2 = current2.min()
on_off_ratio2 = I_on2/I_off2
print(f"On Current I_on ≈ {I_on2:.3e} A ")
print(f"On Current I_off ≈ {I_off2:.3e} A")
print(f"Current on-off ratio ≈ {on_off_ratio2:.2e} ")
I_on3 = current3.max()
I_off3 = current3.min()
on_off_ratio3 = I_on3/I_off3
print(f"On Current I_on ≈ {I_on3:.3e} A ")
print(f"On Current I_off ≈ {I_off3:.3e} A")
print(f"Current on-off ratio ≈ {on_off_ratio3:.2e} ")
I_on4 = current4.max()
I_off4 = current4.min()
on_off_ratio4 = I_on4/I_off4
print(f"On Current I_on ≈ {I_on4:.3e} A ")
print(f"On Current I_off ≈ {I_off4:.3e} A")
print(f"Current on-off ratio ≈ {on_off_ratio4:.2e} ")
I_on5 = current5.max()
I_off5 = current5.min()
on_off_ratio5 = I_on5/I_off5
print(f"On Current I_on ≈ {I_on5:.3e} A ")
print(f"On Current I_off ≈ {I_off5:.3e} A")
print(f"Current on-off ratio ≈ {on_off_ratio5:.2e} ")
I_on6 = current6.max()
I_off6 = current6.min()
on_off_ratio6 = I_on6/I_off6
print(f"On Current I_on ≈ {I_on6:.3e} A ")
print(f"On Current I_off ≈ {I_off6:.3e} A")
print(f"Current on-off ratio ≈ {on_off_ratio6:.2e} ")
I_on7 = current7.max()
I_off7 = current7.min()
on_off_ratio7 = I_on7/I_off7
print(f"On Current I_on ≈ {I_on7:.3e} A ")
print(f"On Current I_off ≈ {I_off7:.3e} A")
print(f"Current on-off ratio ≈ {on_off_ratio7:.2e} ")

# Subthreshold Current (average in low-Vgs region)
subthreshold_current1 = current1[voltage < vth1].mean()
print(f"Subthreshold Current = {subthreshold_current1:.3e} A")
subthreshold_current2 = current2[voltage < vth2].mean()
print(f"Subthreshold Current = {subthreshold_current2:.3e} A")
subthreshold_current3 = current3[voltage < vth3].mean()
print(f"Subthreshold Current = {subthreshold_current3:.3e} A")
subthreshold_current4 = current4[voltage < vth4].mean()
print(f"Subthreshold Current = {subthreshold_current4:.3e} A")
subthreshold_current5 = current5[voltage < vth5].mean()
print(f"Subthreshold Current = {subthreshold_current5:.3e} A")
subthreshold_current6 = current6[voltage < vth6].mean()
print(f"Subthreshold Current = {subthreshold_current6:.3e} A")
subthreshold_current7 = current7[voltage < vth7].mean()
print(f"Subthreshold Current = {subthreshold_current7:.3e} A")

# Effective Mobility
Vds = 0.1
Cox = 6.27e-7
W = 0.5e-4
L = 1e-4
mobility_eff1 = (gm1.max() * L) / (Cox * W * Vds)
print(f"Effective Mobility = {mobility_eff1:.2e} cm²/V·s")
mobility_eff2 = (gm2.max() * L) / (Cox * W * Vds)
print(f"Effective Mobility = {mobility_eff2:.2e} cm²/V·s")
mobility_eff3 = (gm3.max() * L) / (Cox * W * Vds)
print(f"Effective Mobility = {mobility_eff3:.2e} cm²/V·s")
mobility_eff4 = (gm4.max() * L) / (Cox * W * Vds)
print(f"Effective Mobility = {mobility_eff4:.2e} cm²/V·s")
mobility_eff5 = (gm5.max() * L) / (Cox * W * Vds)
print(f"Effective Mobility = {mobility_eff5:.2e} cm²/V·s")
mobility_eff6 = (gm6.max() * L) / (Cox * W * Vds)
print(f"Effective Mobility = {mobility_eff6:.2e} cm²/V·s")
mobility_eff7 = (gm7.max() * L) / (Cox * W * Vds)
print(f"Effective Mobility = {mobility_eff7:.2e} cm²/V·s")


# gm/Id Efficiency
gm_eff1 = gm1 / current1_smooth
gm_eff_mean1 = np.mean(gm_eff1[(voltage > vth1) & (voltage < 10)])
print(f"Mean gm/Id Efficiency = {gm_eff_mean1:.3f} 1/V")
gm_eff2 = gm2 / current2_smooth
gm_eff_mean2 = np.mean(gm_eff2[(voltage > vth2) & (voltage < 10)])
print(f"Mean gm/Id Efficiency = {gm_eff_mean2:.3f} 1/V")
gm_eff3 = gm3 / current3_smooth
gm_eff_mean3 = np.mean(gm_eff3[(voltage > vth3) & (voltage < 10)])
print(f"Mean gm/Id Efficiency = {gm_eff_mean3:.3f} 1/V")
gm_eff4 = gm4 / current4_smooth
gm_eff_mean4 = np.mean(gm_eff4[(voltage > vth4) & (voltage < 10)])
print(f"Mean gm/Id Efficiency = {gm_eff_mean4:.3f} 1/V")
gm_eff5 = gm5 / current5_smooth
gm_eff_mean5 = np.mean(gm_eff5[(voltage > vth5) & (voltage < 10)])
print(f"Mean gm/Id Efficiency = {gm_eff_mean5:.3f} 1/V")
gm_eff6 = gm6 / current6_smooth
gm_eff_mean6 = np.mean(gm_eff6[(voltage > vth6) & (voltage < 10)])
print(f"Mean gm/Id Efficiency = {gm_eff_mean6:.3f} 1/V")
gm_eff7 = gm7 / current7_smooth
gm_eff_mean7 = np.mean(gm_eff7[(voltage > vth7) & (voltage < 10)])
print(f"Mean gm/Id Efficiency = {gm_eff_mean7:.3f} 1/V")

# Turn-on Voltage (peak slope in log scale)
log_ids1 = np.log10(current1_smooth)
log_ids1_smooth = savgol_filter(log_ids1, window_length=155, polyorder=3)
dlog_ids1 = np.gradient(log_ids1_smooth, voltage)
turn_on_index1 = np.argmax(dlog_ids1)
turn_on_voltage1 = voltage.iloc[turn_on_index1]
print(f"Turn-on Voltage = {turn_on_voltage1:.3f} V")

log_ids2 = np.log10(current2_smooth)
log_ids2_smooth = savgol_filter(log_ids2, window_length=145, polyorder=3)
dlog_ids2 = np.gradient(log_ids2_smooth, voltage)
turn_on_index2 = np.argmax(dlog_ids2)
turn_on_voltage2 = voltage.iloc[turn_on_index2]
print(f"Turn-on Voltage = {turn_on_voltage2:.3f} V")

log_ids3 = np.log10(current3_smooth)
log_ids3_smooth = savgol_filter(log_ids3, window_length=135, polyorder=3)
dlog_ids3 = np.gradient(log_ids3_smooth, voltage)
turn_on_index3 = np.argmax(dlog_ids3)
turn_on_voltage3 = voltage.iloc[turn_on_index3]
print(f"Turn-on Voltage = {turn_on_voltage3:.3f} V")

log_ids4 = np.log10(current4_smooth)
log_ids4_smooth = savgol_filter(log_ids4, window_length=125, polyorder=3)
dlog_ids4 = np.gradient(log_ids4_smooth, voltage)
turn_on_index4 = np.argmax(dlog_ids4)
turn_on_voltage4 = voltage.iloc[turn_on_index4]
print(f"Turn-on Voltage = {turn_on_voltage4:.3f} V")

log_ids5 = np.log10(current5_smooth)
log_ids5_smooth = savgol_filter(log_ids5, window_length=115, polyorder=3)
dlog_ids5 = np.gradient(log_ids5_smooth, voltage)
turn_on_index5 = np.argmax(dlog_ids5)
turn_on_voltage5 = voltage.iloc[turn_on_index5]
print(f"Turn-on Voltage = {turn_on_voltage5:.3f} V")

log_ids6 = np.log10(current6_smooth)
log_ids6_smooth = savgol_filter(log_ids6, window_length=105, polyorder=3)
dlog_ids6 = np.gradient(log_ids6_smooth, voltage)
turn_on_index6 = np.argmax(dlog_ids6)
turn_on_voltage6 = voltage.iloc[turn_on_index6]
print(f"Turn-on Voltage = {turn_on_voltage6:.3f} V")

log_ids7 = np.log10(current7_smooth)
log_ids7_smooth = savgol_filter(log_ids7, window_length=105, polyorder=3)
dlog_ids7 = np.gradient(log_ids7_smooth, voltage)
turn_on_index7 = np.argmax(dlog_ids7)
turn_on_voltage7 = voltage.iloc[turn_on_index7]
print(f"Turn-on Voltage = {turn_on_voltage7:.3f} V")