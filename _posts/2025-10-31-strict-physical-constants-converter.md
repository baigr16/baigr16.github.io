---
layout: article
title: 严格物理常数单位换算器
mathjax: true
---

# 严格物理常数单位换算器（含常数与公式显示）

> 输入任意一项后按 **Enter** 或 **空格键**计算。  
> 支持波长 λ (nm)、能量 E (eV / a.u.)、周期 T (fs)、频率 ν (THz)、波数 ν̃ (cm⁻¹) 之间的精确换算。  

---

## 输入与计算

<table>
<tr><th>量</th><th>符号</th><th>输入 / 输出</th></tr>
<tr><td>波长</td><td>λ (nm)</td>
    <td><input id="lambda" onkeypress="if(event.key==='Enter'||event.key===' ') convertFrom('lambda')"></td></tr>
<tr><td>能量</td><td>E (eV)</td>
    <td><input id="energy" onkeypress="if(event.key==='Enter'||event.key===' ') convertFrom('energy')"></td></tr>
<tr><td>周期</td><td>T (fs)</td>
    <td><input id="period" onkeypress="if(event.key==='Enter'||event.key===' ') convertFrom('period')"></td></tr>
<tr><td>频率</td><td>ν (THz)</td>
    <td><input id="freq" onkeypress="if(event.key==='Enter'||event.key===' ') convertFrom('freq')"></td></tr>
<tr><td>波数</td><td>ν̃ (cm⁻¹)</td>
    <td><input id="wavenum" onkeypress="if(event.key==='Enter'||event.key===' ') convertFrom('wavenum')"></td></tr>
<tr><td>能量</td><td>E (a.u.)</td>
    <td><input id="eau" onkeypress="if(event.key==='Enter'||event.key===' ') convertFrom('eau')"></td></tr>
</table>

---

## 常数与换算

<pre id="constants"></pre>

---

## 换算公式

<pre>
E = hc / λ
T = h / E
ν = 1 / T
ν̃ = 1 / λ
E(a.u.) = E(eV) / 27.211386
1 eV = 8065.544 cm⁻¹ = 241.799 THz = 0.0367493 a.u.
hc = 1239.841984 eV·nm
</pre>

---

<script>
// ======== 常数 ========
const _PLANCK   = 0.6626075540000000E-33;     // h [J·s]
const _CLIGHT   = 29979245800.00000;         // c [cm/s]
const _JOULE    = 0.4359748200000000E-17;    // J per Hartree
const _TOEV     = 27.21138505000000;         // eV per Hartree
const _TOHERTZ  = 6579683920729000.0;        // Hz per Hartree
const _EV_HARTREE = 0.3674932379085202E-01;  // Hartree per eV

// ======== 派生常数 ========
const hc_J_cm = _PLANCK * _CLIGHT;
const hc_eV_nm = hc_J_cm * (_TOEV / _JOULE) * 1e7;
const eV_to_Hz = _TOHERTZ * _EV_HARTREE;
const eV_to_THz = eV_to_Hz / 1e12;
const J_per_eV = _JOULE / _TOEV;
const h_eV_fs = _PLANCK / J_per_eV * 1e15;
const eV_to_cm1 = 1e7 / hc_eV_nm;

// ======== 格式化输出 ========
function fmt(x,n=6){return (isFinite(x)? x.toExponential(n): '');}

// ======== 转换函数 ========
function convertFrom(type){
  const v = parseFloat(document.getElementById(type).value);
  if (!isFinite(v)) return;

  let lambda_nm, E_eV, T_fs, nu_THz, wavenum_cm1, E_au;

  if (type === 'lambda'){
    lambda_nm = v;
    E_eV = hc_eV_nm / lambda_nm;
  } else if (type === 'energy'){
    E_eV = v;
    lambda_nm = hc_eV_nm / E_eV;
  } else if (type === 'period'){
    T_fs = v;
    E_eV = h_eV_fs / T_fs;
    lambda_nm = hc_eV_nm / E_eV;
  } else if (type === 'freq'){
    nu_THz = v;
    E_eV = (nu_THz * 1e12) / eV_to_Hz;
    lambda_nm = hc_eV_nm / E_eV;
  } else if (type === 'wavenum'){
    wavenum_cm1 = v;
    E_eV = wavenum_cm1 / eV_to_cm1;
    lambda_nm = hc_eV_nm / E_eV;
  } else if (type === 'eau'){
    E_au = v;
    E_eV = E_au * _TOEV;
    lambda_nm = hc_eV_nm / E_eV;
  }

  // 衍生计算
  nu_THz = E_eV * eV_to_THz;
  T_fs = h_eV_fs / E_eV;
  wavenum_cm1 = E_eV * eV_to_cm1;
  E_au = E_eV / _TOEV;

  // 输出更新
  document.getElementById('lambda').value  = lambda_nm.toFixed(6);
  document.getElementById('energy').value  = E_eV.toFixed(9);
  document.getElementById('period').value  = T_fs.toFixed(6);
  document.getElementById('freq').value    = nu_THz.toFixed(6);
  document.getElementById('wavenum').value = wavenum_cm1.toFixed(6);
  document.getElementById('eau').value     = E_au.toFixed(9);
}

// ======== 显示常数 ========
document.getElementById("constants").innerText =
`h  = ${fmt(_PLANCK)} J·s
c  = ${fmt(_CLIGHT)} cm/s
1 Hartree = ${_TOEV} eV = ${(1/_EV_HARTREE).toFixed(6)} eV
1 eV = ${fmt(eV_to_cm1,3)} cm⁻¹ = ${fmt(eV_to_THz,3)} THz = ${(1/_TOEV).toFixed(8)} a.u.
hc = ${hc_eV_nm.toFixed(6)} eV·nm
h  = ${(h_eV_fs).toFixed(6)} eV·fs`;
</script>
