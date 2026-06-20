#--------PART A---------
# אינטגרציה והתאמה לינארית

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import re

# === פרמטרים פיזיקליים לחישוב ===
I0 = 1.0       # עוצמת מקור ה-LED (יש לעדכן מהאינטגרציה על ספקטרום המקור)
l_path = 1.0   # אורך המסלול האופטי בקיווטה [cm]
delta_lambda = 0.0 # שגיאת הרזולוציה/אורך הגל של הספקטרומטר [nm] 

quantum_yields = {
    "Fluorescein": 0.97,
    "Rhodamine 6G": 0.95,
    "Rhodamine B": 0.70
}

materials_files = {
    "Rhodamine 6G": ["6G_0.0025mm.ods", "6G_0.005mm.ods", "6G_0.01mm.ods", "6G_0.025mm.ods", "6G_0.05mm.ods", "6G_0.1mm.ods"],
    "Rhodamine B": ["B_0.0025mm.ods", "B_0.005MM.ods", "B_0.01mm.ods", "B_0.025mm.ods", "B_0.05mm.ods", "RHBB_0.1mm.ods"],
    "Fluorescein": ["fe_0.0025mm.ods", "fe_0.005mm.ods", "fe_0.01mm.ods", "fe_0.025mm.ods", "FE_0.05mm.ods", "fe_0.1mm.ods"]
}

min_wave = 480
max_wave = 650

def linear_model(x, m, b):
    return m * x + b

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

for idx, (material_name, files) in enumerate(materials_files.items()):
    concentrations = []
    integrals = []
    integrals_errors = [] 
    
    for f in files:
        try:
            df = pd.read_excel(f, engine="odf")
            
            if len(df.columns) == 1 and ';' in df.columns[0]:
                col_name = df.columns[0]
                new_cols = col_name.split(';')
                df = df[col_name].str.split(';', expand=True)
                df.columns = new_cols

            for col in df.columns:
                df[col] = pd.to_numeric(df[col], errors='coerce')

            x_col = df.columns[0]
            y_col = df.columns[1]
            df = df.dropna(subset=[x_col, y_col])

            mask = (df[x_col] >= min_wave) & (df[x_col] <= max_wave)
            df_filtered = df[mask]
            
            if len(df_filtered) > 1:
                x_val = df_filtered[x_col].values
                y_val = df_filtered[y_col].values
                
                # אינטגרציה נומרית (כלל הטרפז)
                dx = np.diff(x_val)
                y_mid = (y_val[:-1] + y_val[1:]) / 2.0
                area = np.sum(y_mid * dx)
                
                # חישוב שגיאות (Shot Noise + Wavelength Error)
                sigma_y = np.sqrt(np.abs(y_val))
                sigma_y_mid = np.sqrt(sigma_y[:-1]**2 + sigma_y[1:]**2) / 2.0
                sigma_dx = np.sqrt(2) * delta_lambda
                sigma_area_i = np.sqrt((dx * sigma_y_mid)**2 + (y_mid * sigma_dx)**2)
                err_area = np.sqrt(np.sum(sigma_area_i**2))

                match = re.search(r'(\d*\.?\d+)mm', f, re.IGNORECASE)
                if match:
                    c_val = np.round(float(match.group(1)), 4)
                    concentrations.append(c_val)
                    integrals.append(area)
                    integrals_errors.append(err_area)
                    
        except Exception as e:
            print(f"Error processing {f}: {e}")

    if concentrations:
        sorted_indices = np.argsort(concentrations)
        c_arr = np.array(concentrations)[sorted_indices]
        i_arr = np.array(integrals)[sorted_indices]
        i_err_arr = np.array(integrals_errors)[sorted_indices]
        
        color = colors[idx % len(colors)]
        
        # ---------------------------------------------------------
        # גרף 1: כל טווח הריכוזים עם צלבי שגיאה
        # ---------------------------------------------------------
        # fmt='-o' יוצר קו רציף בין הנקודות יחד עם סמנים עגולים
        ax1.errorbar(c_arr, i_arr, yerr=i_err_arr, fmt='-o', color=color, alpha=0.6, capsize=4, label=material_name)
        
        # ---------------------------------------------------------
        # גרף 2: התאמה לינארית ל-3 הנקודות הראשונות עם סורגי שגיאה
        # ---------------------------------------------------------
        if len(c_arr) >= 3:
            c_fit = c_arr[:3]
            i_fit = i_arr[:3]
            i_err_fit = i_err_arr[:3]
            
            ax2.errorbar(c_fit, i_fit, yerr=i_err_fit, fmt='o', color=color, capsize=4, label=f"{material_name} (Data)")
            
            popt, pcov = curve_fit(linear_model, c_fit, i_fit, sigma=i_err_fit, absolute_sigma=False)
            slope, intercept = popt
            slope_err = np.sqrt(np.diag(pcov))[0] 
            
            fit_line = slope * c_fit + intercept
            ax2.plot(c_fit, fit_line, linestyle='--', color=color, linewidth=2, label=f"{material_name} Fit")
            
            phi_f = quantum_yields.get(material_name, 1.0)
            scaling_factor = (2.303 * I0 * phi_f * l_path)
            
            k_epsilon = slope / scaling_factor
            k_epsilon_err = slope_err / scaling_factor
            
            print(f"--- {material_name} ---")
            print(f"Linear Slope (m): {slope:.4f} ± {slope_err:.4f}")
            print(f"k * Epsilon: {k_epsilon:.4f} ± {k_epsilon_err:.4f} cm^-1 mM^-1\n")

# עיצוב הגרף הראשון (כללי)
ax1.set_xlabel("Concentration [mM]", fontsize=12)
ax1.set_ylabel("Integrated Intensity [a.u.]", fontsize=12)
ax1.set_title("Integrated Intensity vs Concentration", fontsize=13)
ax1.legend(fontsize=10)
ax1.grid(True, linestyle='--', alpha=0.7)

# עיצוב הגרף השני (3 נקודות ראשונות)
ax2.set_xlabel("Concentration [mM]", fontsize=12)
ax2.set_ylabel("Integrated Intensity [a.u.]", fontsize=12)
ax2.set_title("Linear Fit of the First 3 Concentrations (with Errors)", fontsize=13)
ax2.legend(fontsize=10)
ax2.grid(True, linestyle='--', alpha=0.7)
ax2.axhline(0, color='black', linewidth=0.8, linestyle=':')

plt.tight_layout()
plt.show()

# Data plot 
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid
import re

# רשימת הקבצים באותה התיקייה
#files = [
 #   "B_0.0025mm.ods",
  #  "B_0.005MM.ods",
   # "B_0.01mm.ods",
    #"B_0.025mm.ods",
    #"B_0.05mm.ods",
    #"RHBB_0.1mm.ods"
#]

files = [
    "6G_0.0025mm.ods",
    "6G_0.005mm.ods",
    "6G_0.01mm.ods",
    "6G_0.025mm.ods",
    "6G_0.05mm.ods",
    "6G_0.1mm.ods"
]
# טווח האינטגרציה לפלורסנציה (ניתן לשנות כאן כדי להוציא את המנורה)
# נניח שהמנורה מסתיימת ב-470nm בערך, נבצע אינטגרציה מ-480 עד 800.
min_wave = 480
max_wave = 650

thicknesses = []
integrals = []

plt.figure(figsize=(10, 6))

for f in files:
    # קריאת הקובץ (שים לב: נדרש להתקין odfpy בעזרת pip install odfpy)
    df = pd.read_excel(f, engine="odf")

    # פיצול עמודות במידה והנתונים מחוברים בנקודה-פסיק כפי שמופיע בקבצים שלך
    if len(df.columns) == 1 and ';' in df.columns[0]:
        col_name = df.columns[0]
        new_cols = col_name.split(';')
        df = df[col_name].str.split(';', expand=True)
        df.columns = new_cols

    # המרת העמודות למספרים (התעלמות משגיאות במקרה של שורות טקסט)
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce')

    # הנחה שהעמודה הראשונה היא אורך הגל (X) והשנייה עוצמה (Y)
    x_col = df.columns[0]
    y_col = df.columns[1]
    df = df.dropna(subset=[x_col, y_col])

    # ציור הגרף (עוצמה כתלות באורך גל לכל דגימה)
    plt.plot(df[x_col], df[y_col], label=f)

    # יצירת המסיכה לאינטגרציה בטווח הרצוי בלבד
    mask = (df[x_col] >= min_wave) & (df[x_col] <= max_wave)
    df_filtered = df[mask]

    # ביצוע האינטגרציה (Trapezoidal Rule)
    area = np.trapz(df_filtered[y_col].values, df_filtered[x_col].values)

    # חילוץ עובי הדגימה מתוך שם הקובץ (למשל שלפת המספר 0.05 מתוך FE_0.05)
    match = re.search(r'_(\d*\.?\d+)mm', f, re.IGNORECASE)
    if match:
      concentration = float(match.group(1))
      thicknesses.append(concentration)
      integrals.append(area)

# עיצוב גרף הספקטרומים הראשון
plt.xlim(400, 740) # ניתן להרחיב לכל הטווח
plt.axvline(min_wave, color='gray', linestyle='--', label='Integration Start')
plt.axvline(max_wave, color='gray', linestyle=':', label='Integration End')
plt.xlabel("Wavelength [nm]")
plt.ylabel("Intensity [a.u.]")
plt.title("Intensity vs Wavelength")
plt.legend()
plt.show()

# סידור נתוני העובי מהקטן לגדול כדי שהקו בגרף יצא רציף
sorted_indices = np.argsort(thicknesses)
thicknesses = np.array(thicknesses)[sorted_indices]
integrals = np.array(integrals)[sorted_indices]

# יצירת הגרף השני: עוצמה שעברה אינטגרציה (שטח הפלורסנציה) כתלות בעובי
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# חלק ראשון - ציר Y ליניארי
axes[0].plot(thicknesses, integrals, marker='o', linestyle='-', color='b')
axes[0].set_xlabel("Concentration [mM]")
axes[0].set_ylabel("Integrated Fluorescence Intensity [a.u.]")
axes[0].set_title("Integrated Intensity vs Concentration (Linear)")
axes[0].grid(True)

# חלק שני - ציר Y לוגריתמי
axes[1].plot(thicknesses, integrals, marker='o', linestyle='-', color='r')
axes[1].set_yscale("log")
axes[1].set_xlabel("Concentration [mM]")
axes[1].set_ylabel("Integrated Fluorescence Intensity (Log scale) [a.u.]")
axes[1].set_title("Integrated Intensity vs Concentration (Log Y)")
axes[1].grid(True, which="both", ls="--")

plt.show()


#--------PART B---------
# המרת תמונה לפלוט של עוצמה כתלות באורך
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

filename = 'RHBB 0.0025mM.jpeg'

# קריאת התמונה
try:
    A = mpimg.imread(filename)
except FileNotFoundError:
    print("Please make sure to upload the image to the Colab environment.")
    A = np.zeros((1536, 2048, 3))

# המרה ל-double (ערכים בין 0 ל-1)
if A.dtype == np.uint8:
    A1 = A.astype(float) / 255.0
else:
    A1 = A.copy()

height, width = A1.shape[0], A1.shape[1]
print(f"Image dimensions: Height = {height}, Width = {width}")

# יצירת תמונה המהווה ממוצע של כלל הערוצים (RGB)
A_mean = np.mean(A1, axis=2)

# =======================================================
#           בחירת קואורדינטות (חצי אוטומטי)
# =======================================================

# 1. זיהוי אוטומטי של ציר Y (השורה המרכזית של הקרן המוארת)
row_profiles = np.mean(A_mean, axis=1)
row_index = np.argmax(row_profiles)

# 2. הגדרה ידנית של ציר X (גבולות תחילת וסיום התמיסה)
start_col = 600
end_col = 1500

# הגנה מפני חריגה מגבולות התמונה
start_col = min(start_col, width - 1)
end_col = min(end_col, width)

print(f"\n[AUTO-DETECTED] Row index (Y): {row_index}")
print(f"[MANUAL] Columns (X): From {start_col} to {end_col}\n")

# ---- חלון ראשון: התמונה עם קו הסימון ----
plt.figure("image", figsize=(10, 5))
plt.imshow(A_mean, cmap='viridis', aspect='auto')
plt.plot([start_col, end_col], [row_index, row_index], color='red', linewidth=2, label='Selection Line')
plt.colorbar(label='Intensity')
plt.title(f'{filename}')
plt.legend()

# שמירת התמונה הראשונה
plt.savefig(f'image_selection {filename[:-5]}.png', bbox_inches='tight')
print("Saved image selection to 'image_selection.png'")
plt.show()

# ---- חילוץ ועיבוד הנתונים לגרף ----
Av = np.mean(A1[row_index, start_col:end_col, :], axis=1)
Av = Av[::-1] # הפיכת הכיוון לקבלת אקספוננט דועך

x = np.arange(0, len(Av))

# חישוב הלוגריתם של העוצמה
with np.errstate(divide='ignore', invalid='ignore'):
    Avl = np.log(Av)

# סינון ערכים לא תקינים
valid = np.isfinite(x) & np.isfinite(Avl)
x_fit = x[valid]
Avl_fit = Avl[valid]

# ---- ביצוע התאמה ליניארית (Linear Fit) ----
p, cov = np.polyfit(x_fit, Avl_fit, 1, cov=True)
a, b = p
perr = np.sqrt(np.diag(cov))
a_err, b_err = perr

# יצירת מחרוזת הטקסט עם התוצאות
results_text = f"""========================================
             FIT RESULTS
========================================
Formula: ln(Intensity) = a * x + b
Slope (a)       = {a:.4e} ± {a_err:.4e}
Intercept (b)   = {b:.4e} ± {b_err:.4e}
========================================"""

# הדפסה למסך
print(results_text)

# שמירת התוצאות לקובץ טקסט
with open(f'fit_results {filename[:-5]}.txt', 'w', encoding='utf-8') as f:
    f.write(results_text)
print("Saved fit results to 'fit_results.txt'")

# חישוב קו ההתאמה והשארים
fit_line = a * x + b
residuals = Avl - fit_line

# ---- חלון שני: גרף העוצמה, ההתאמה והשארים ----
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True,
                               gridspec_kw={'height_ratios': [3, 1]})

# גרף עליון
ax1.plot(x, Avl, 'o', markersize='2', label='Logarithmic Intensity (Data)', color='blue')
ax1.plot(x, fit_line, '-', label='Linear Fit', color='red', linewidth=2)
ax1.set_ylabel('log(Intensity)')
ax1.set_title('Linear Fit of Intensity vs Position')
ax1.grid(True, which='both', linestyle='--', alpha=0.5)
ax1.legend()

# גרף תחתון
ax2.plot(x, residuals, 'o', markersize='2', color='blue', label='Residuals ($\Delta y$)')
ax2.axhline(0, color='black', linestyle='--', linewidth=1)
ax2.set_xlabel('Position x [pixels]')
ax2.set_ylabel('Residuals')
ax2.grid(True, linestyle='--', alpha=0.5)
ax2.legend()

plt.tight_layout()

# שמירת התמונה השנייה
plt.savefig(f'fit_and_residuals {filename[:-5]}.png', bbox_inches='tight')
print("Saved plots to 'fit_and_residuals.png'")
plt.show()



# מייצר גרף מאוחד של כל הריכוזים
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import re

# =======================================================
# 1. הגדרות כיול מרחבי (Spatial Calibration)
# =======================================================
PIXELS_PER_CM = 150

# =======================================================
# 2. נתוני התמונות (קלט המשתמש)
# =======================================================
image_data = [
    {"filename": "RHBB 0.0025mM.jpeg", "start_col": 600, "end_col": 1500},
    {"filename": "RHBB 0.005mM.jpeg",  "start_col": 600, "end_col": 1500},
    {"filename": "RHBB 0.01mM.jpeg",   "start_col": 500, "end_col": 1550},
    {"filename": "RHBB 0.025mM.jpeg",  "start_col": 500, "end_col": 1550},
    {"filename": "RHBB 0.05mM.jpeg",   "start_col": 500, "end_col": 1550},
    {"filename": "RHBB 0.1mM.jpeg",    "start_col": 500, "end_col": 1550}
]

# רשימה לשמירת תוצאות הטבלה
results_table = []

# פתיחת גרף משותף להצגת כלל התוצאות
plt.figure(figsize=(12, 8))

for data in image_data:
    filename = data["filename"]
    start_col = data["start_col"]
    end_col = data["end_col"]

    # קריאת התמונה
    try:
        A = mpimg.imread(filename)
    except FileNotFoundError:
        print(f"[Warning] File not found: {filename}. Skipping...")
        continue

    # המרה לערכים צפים (0-1)
    if A.dtype == np.uint8:
        A1 = A.astype(float) / 255.0
    else:
        A1 = A.copy()

    height, width = A1.shape[:2]

    # זיהוי אוטומטי של ציר Y (מישור הלייזר/הקרן המוארת)
    A_mean = np.mean(A1, axis=2)
    row_profiles = np.mean(A_mean, axis=1)
    row_index = np.argmax(row_profiles)

    # הגנה מפני חריגה מגבולות
    start_col = max(0, min(start_col, width - 1))
    end_col = max(0, min(end_col, width))

    # חילוץ חתך העוצמה והפיכת כיוון (כך שהקרן דועכת משמאל לימין)
    Av = np.mean(A1[row_index, start_col:end_col, :], axis=1)
    Av = Av[::-1]

    # יצירת ציר X והמרה לסנטימטרים
    x_pixels = np.arange(len(Av))
    x_cm = x_pixels / PIXELS_PER_CM

    # =======================================================
    # נרמול לפי העוצמה המקסימלית ומעבר ללוגריתם
    # =======================================================
    # חלוקה בעוצמה המקסימלית כך שהמקסימום יהיה 1
    Av_norm_linear = Av / np.max(Av)

    # מעבר ללוגריתם טבעי (הערך המקסימלי יהיה עכשיו בדיוק 0)
    with np.errstate(divide='ignore', invalid='ignore'):
        Avl_norm = np.log(Av_norm_linear)

    # סינון ערכים לא תקינים (NaN / Inf) לפני ההתאמה
    valid = np.isfinite(x_cm) & np.isfinite(Avl_norm)
    x_fit = x_cm[valid]
    Avl_fit = Avl_norm[valid]

    if len(x_fit) == 0:
        print(f"Skipping {filename}: Not enough valid data points.")
        continue

    # ביצוע התאמה ליניארית לנתונים המנורמלים: y = ax + b
    p, cov = np.polyfit(x_fit, Avl_fit, 1, cov=True)
    a, b = p

    # חישוב אלפא (השיפוע השלילי)
    alpha = -a

    # חילוץ הריכוז מהשם של הקובץ, חישוב אפסילון והגדרת שם במקרא
    match = re.search(r'([0-9.]+)\s*mM', filename)
    if match:
        conc_mM = float(match.group(1))
        conc_M = conc_mM * 1e-3 # המרה למולר (M)

        # חישוב epsilon [M^-1 cm^-1]
        epsilon = alpha / (conc_M * np.log(10))
        results_table.append([filename, conc_mM, alpha, epsilon])

        # הגדרת התווית שתופיע במקרא של הגרף (רק הריכוז)
        label_str = f"{match.group(1)} mM"
    else:
        results_table.append([filename, None, alpha, None])
        # גיבוי במקרה שהקובץ לא מכיל ריכוז בפורמט המוכר
        label_str = filename.split('.')[0]

    # יצירת קו ההתאמה
    fit_line_norm = a * x_fit + b

    # ציור הנקודות המנורמלות ושמירת הצבע
    scatter_plot = plt.plot(x_fit, Avl_fit, 'o', markersize=2, alpha=0.3)
    color = scatter_plot[0].get_color()

    # הצגת קו ההתאמה הליניארית באותו צבע עם התווית החדשה
    plt.plot(x_fit, fit_line_norm, '-', linewidth=2, color=color, label=label_str)

# =======================================================
# 3. הדפסת טבלת התוצאות לקונסולה
# =======================================================
print("="*80)
print(f"{'Filename':<25} | {'Conc [mM]':<10} | {'alpha [cm^-1]':<15} | {'epsilon [M^-1 cm^-1]':<20}")
print("-" * 80)
for row in results_table:
    if row[1] is not None:
        print(f"{row[0]:<25} | {row[1]:<10.4f} | {row[2]:<15.4f} | {row[3]:<20.2f}")
    else:
        print(f"{row[0]:<25} | {'N/A':<10} | {row[2]:<15.4f} | {'N/A':<20}")
print("="*80)

# עיצוב הגרף המשותף
plt.xlabel('Position $x$ [cm]')
plt.ylabel(r'Normalized $\ln(I / I_{max})$')
plt.title('Normalized intensity as a function of position - Rhodamine B')
plt.grid(True, which='both', linestyle='--', alpha=0.6)

# הזזת מקרא התמונות מחוץ לגרף
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.show()