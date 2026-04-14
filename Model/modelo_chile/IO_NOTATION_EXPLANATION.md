# Input-Output Matrix Notation Explanation

## The Confusion

You're right to be confused! The model uses **different index ordering** in different equations. Here's why:

---

## Convention Used

Let's define:
- **i** = producing/supplying sector  
- **j** = using/purchasing sector

---

## 1. DATA FORMAT: IO_2021_chile.csv

**Format: (i,j) where i = supplier, j = purchaser**

```
Data(i,j) = amount that sector j purchases from sector i
```

Example from CSV:
- Row 1 (Agriculture), Column 3 (Manufacturing) = 0.001192
- This means: **Manufacturing (j=3) purchases 0.001192 from Agriculture (i=1)**

After normalization by columns:
```
betax(i,j) = Data(i,j) / sum_i Data(i,j)  
          = share that sector j purchases from sector i
```

---

## 2. MODEL PARAMETER: modbeta in main_SOE.m

**After transpose: modbeta(j,i) where i = supplier, j = purchaser**

```matlab
modbeta = betax'  % Transpose
```

So: `modbeta(j,i) = betax(i,j)`

This means:
```
modbeta(j,i) = share that sector j purchases from sector i
```

**Using your convention (i=producer, j=user):**
- modbeta has indices in order: (purchaser, supplier)
- modbeta has indices in order: (j, i)  ✓

---

## 3. MODEL CODE: NK_IOSOE_lev.mod

The .mod file assigns: `beta_@{i}_@{j} = modbeta(@{i},@{j})`

### Equation 1: Price Index (Line 291-295)

```
PM_i = [ sum_j beta_{i,j} * (PH_j)^(1-eps) ]^(1/(1-eps))
```

**Interpretation:**  
- PM_i = price index for intermediate inputs used by sector i
- beta_{i,j} is paired with PH_j (price of sector j's output)
- This means: **sector i is PURCHASING from sector j**

**Using your notation (i=producer, j=user):**  
Here we need i=purchaser, j=supplier, so the indices are **(j,i)** in your convention!

But the code writes it as `beta_{i,j}`, which means:
- First index (i) = **purchaser** 
- Second index (j) = **supplier**
- This is opposite to your convention!

### Equation 2: Market Clearing (Line 311-313)

```
Y_i = consumption + exports + sum_j beta_{j,i} * M_j
```

**Interpretation:**  
- Y_i = output produced by sector i
- beta_{j,i} weights how much sector j uses from sector i
- This means: **sector j is PURCHASING from sector i**

**Using your notation (i=producer, j=user):**  
Here i=supplier, j=purchaser, so indices are **(j,i)** = (user, producer) ✓  
This matches your convention!

---

## SUMMARY: Why Different Orderings?

The model uses **TWO DIFFERENT INTERPRETATIONS** depending on the equation:

| Location | Code Notation | Means | Convention (i=prod, j=use) |
|----------|---------------|-------|----------------------------|
| **Price Index** | beta_{i,j} | i purchases from j | (j,i) - REVERSED! |
| **Market Clearing** | beta_{j,i} | j purchases from i | (j,i) - MATCHES! |
| **Data Table** | Table(i,j) | j purchases from i | (j,i) - MATCHES! |
| **modbeta Matrix** | modbeta(a,b) | a purchases from b | - |

---

## The Key Insight

The model parameter **beta_{i,j}** in the .mod file should be read as:

> **"beta with first index i and second index j means i PURCHASES from j"**

This is the **OPPOSITE** of the natural IO table convention where (i,j) means i supplies to j!

That's why we need the transpose in main_SOE.m:
- Data has: (supplier, purchaser)  
- modbeta has: (purchaser, supplier) ← REVERSED for the model equations

---

## Verification

From NK_IOSOE_lev.mod line 134:
```
beta_@{i}_@{j} = modbeta(@{i},@{j})
```

So: `beta_{1,3} = modbeta(1,3)`

From our processing:
- modbeta(1,3) = betax(3,1) = Data(3,1) / sum_i Data(i,1)
- Data(3,1) = row 3 (Manufacturing), col 1 (Agriculture) = 0.006498

This means:
- beta_{1,3} = share that Agriculture (sector 1) purchases from Manufacturing (sector 3)
- Value = 0.006498 ✓

Check in model equation PM_1 (price index for Agriculture):
```
PM_1 = [ ... + beta_{1,3} * (PH_3)^(1-eps) + ... ]^(1/(1-eps))
```
This says Agriculture's input price depends on Manufacturing's price (PH_3), weighted by beta_{1,3}.
This confirms: **beta_{1,3} = how much Agriculture buys from Manufacturing** ✓

---

## Conclusion

**You are correct!** The indices appear in different orders:

1. **Data table**: (i,j) = (producer, user) - natural convention
2. **Model beta_{i,j}**: (i,j) = (user, producer) - reversed!
3. **Model beta_{j,i}** in market clearing: (j,i) = (user, producer) - natural!

The transpose in main_SOE.m exists specifically to convert from the natural data convention to the reversed convention that the price index equation requires.
