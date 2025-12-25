# Graph.py
import numpy as np
import matplotlib.pyplot as plt

# Чтение данных
data_h = np.loadtxt('harmonic_signal.txt', skiprows=1)
idx_h = data_h[:, 0]
orig_h = data_h[:, 1]
filt_h = data_h[:, 2]

data_pw = np.loadtxt('piecewise_signal.txt', skiprows=1)
idx_pw = data_pw[:, 0]
orig_pw = data_pw[:, 1]
filt_pw = data_pw[:, 2]

# Настройка графиков
plt.rcParams['figure.figsize'] = [10, 4]
plt.rcParams['lines.linewidth'] = 1.2
plt.rcParams['axes.grid'] = True
plt.rcParams['grid.alpha'] = 0.3
plt.rcParams['grid.linestyle'] = '--'

# График 1: Гармонический сигнал
plt.figure(figsize=(12, 5))
plt.plot(idx_h[:500], orig_h[:500], '#1f77b4', label='Исходный', alpha=0.8, linewidth=1.5)
plt.plot(idx_h[:500], filt_h[:500], '#ff7f0e', label='Отфильтрованный', alpha=0.9, linewidth=1.5)
plt.xlabel('Отсчеты, j', fontsize=12)
plt.ylabel('Амплитуда', fontsize=12)
plt.title('Гармонический сигнал: DFT ', fontsize=14, fontweight='bold', pad=15)
plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()
plt.savefig('harmonic.png', dpi=150, bbox_inches='tight')
plt.show()

# График 2: Кусочно-постоянный сигнал
plt.figure(figsize=(12, 5))
plt.plot(idx_pw[:500], orig_pw[:500], '#1f77b4', label='Исходный', alpha=0.8, linewidth=1.5)
plt.plot(idx_pw[:500], filt_pw[:500], '#ff7f0e', label='Отфильтрованный', alpha=0.9, linewidth=1.5)
plt.xlabel('Отсчеты, j', fontsize=12)
plt.ylabel('Амплитуда', fontsize=12)
plt.title('Кусочно-постоянный сигнал: DFT ', fontsize=14, fontweight='bold', pad=15)
plt.legend(loc='upper right', fontsize=11)
plt.tight_layout()
plt.savefig('piecewise.png', dpi=150, bbox_inches='tight')
plt.show()

# Дополнительно: Раздельные графики (как на изображениях)
fig, axes = plt.subplots(2, 2, figsize=(14, 8))

# Гармонический сигнал - исходный
axes[0, 0].plot(idx_h[:500], orig_h[:500], '#1f77b4', linewidth=1.5)
axes[0, 0].set_xlabel('Отсчеты, j')
axes[0, 0].set_ylabel('Амплитуда')
axes[0, 0].set_title('Гармонический сигнал: DFT \nИсходный', fontweight='bold')
axes[0, 0].grid(True, alpha=0.3, linestyle='--')
axes[0, 0].set_ylim(-3, 3)

# Гармонический сигнал - отфильтрованный
axes[0, 1].plot(idx_h[:500], filt_h[:500], '#ff7f0e', linewidth=1.5)
axes[0, 1].set_xlabel('Отсчеты, j')
axes[0, 1].set_ylabel('Амплитуда')
axes[0, 1].set_title('Гармонический сигнал: DFT \nОтфильтрованный', fontweight='bold')
axes[0, 1].grid(True, alpha=0.3, linestyle='--')
axes[0, 1].set_ylim(-3, 3)

# Кусочно-постоянный сигнал - исходный
axes[1, 0].plot(idx_pw[:500], orig_pw[:500], '#1f77b4', linewidth=1.5)
axes[1, 0].set_xlabel('Отсчеты, j')
axes[1, 0].set_ylabel('Амплитуда')
axes[1, 0].set_title('Кусочно-постоянный сигнал: DFT \nИсходный', fontweight='bold')
axes[1, 0].grid(True, alpha=0.3, linestyle='--')
axes[1, 0].set_ylim(-0.5, 3)

# Кусочно-постоянный сигнал - отфильтрованный
axes[1, 1].plot(idx_pw[:500], filt_pw[:500], '#ff7f0e', linewidth=1.5)
axes[1, 1].set_xlabel('Отсчеты, j')
axes[1, 1].set_ylabel('Амплитуда')
axes[1, 1].set_title('Кусочно-постоянный сигнал: DFT \nОтфильтрованный', fontweight='bold')
axes[1, 1].grid(True, alpha=0.3, linestyle='--')
axes[1, 1].set_ylim(-0.5, 3)

plt.tight_layout()
plt.savefig('comparison_separated.png', dpi=150, bbox_inches='tight')
plt.show()

print("Графики сохранены:")
print("1. harmonic.png - гармонический сигнал")
print("2. piecewise.png - кусочно-постоянный сигнал")
print("3. comparison_separated.png - раздельные графики как в примере")