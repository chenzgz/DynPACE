
# DynPACE

An R package for dynamic prediction in survival analysis, implementing a framework that combines landmark analysis with Principal Analysis by Conditional Expectation (PACE). Provides both dynamic Cox and RMST modeling approaches and generates personalized dynamic risk prediction curves for individual patients.

## Installation

You can install the development version of DynPACE from [GitHub](https://github.com/chengzhengzhi/DynPACE) with:

```r
# install.packages("devtools")
devtools::install_github("chengzhengzhi/DynPACE")
```

## Main Functions

- `Cox_pace()` - Dynamic Cox prediction model with PACE
- `RMST_pace()` - Dynamic RMST prediction model with PACE  
- `Cox_risk_curve()` - Individual dynamic risk curves for Cox model
- `RMST_risk_curve()` - Individual dynamic risk curves for RMST model

## Data

- `pbc2_example` - Processed subset of Primary Biliary Cholangitis data

## License

GPL-3
