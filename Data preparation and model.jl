#Scenario 02: Limitating the Gas
# - set a CO2 price as c_prod for SMR to 0.0025$, 0.00325$ and 0.004$. (30€ for 1kg CO2)
# - [set the gas capacity from 10^17 to 10^7 MJ] DIDN'T DO
# Funktioniert mit verschiedenen Längen für timesteps
# Regionen haben eigenen Speicher, eigene Produktion, eigene Nachfrage
# Es kann zwischen Regionen zu festlegbaren Kosten transportiert werden
# Alle Regionen greifen auf die gleichen Ressourcenpreise zu
# Es ab und zu wenig produziert, sodass die Nachfrage nicht gedeckt wäre, aber das ist im e^-11 oder so-Bereich

x=1
using Pkg
Pkg.add("Cbc")
Pkg.add("JuMP")
Pkg.add("CSV")
Pkg.add("DataFrames")
Pkg.add("Plots")
Pkg.add("Plotly")
Pkg.add("PlotlyJS")
Pkg.add("StatsPlots")
Pkg.add("PyPlot")
using Cbc, JuMP, CSV, DataFrames, Plots, StatsPlots, PlotlyJS


#read dataframes
df_AnnualPowerDemand = CSV.File(joinpath(@__DIR__, "AnnualPowerDemand.csv")) |> DataFrame
df_AnnualHydrogenDemand = CSV.File(joinpath(@__DIR__, "AnnualHydrogenDemand.csv"), delim=";") |> DataFrame
df_PowerDemandProfile = CSV.File(joinpath(@__DIR__, "PowerDemandProfile.csv")) |> DataFrame
df_WindOnshoreTimeseries = CSV.File(joinpath(@__DIR__, "WindOnshoreTimeseries.csv")) |> DataFrame
df_SolarTimeseries = CSV.File(joinpath(@__DIR__, "SolarTimeseries.csv")) |> DataFrame
df_WindOffshoreTimeseries = CSV.File(joinpath(@__DIR__, "WindOffshoreTimeseries.csv"), type=String, types=Dict(1 => Int)) |> DataFrame
df_CapacityInstalledGW = CSV.File(joinpath(@__DIR__, "capacity_installed_gw.csv")) |> DataFrame
df_AnnualHydrogenDemandRegions = CSV.File(joinpath(@__DIR__, "AnnualHydrogenDemandRegions.csv")) |> DataFrame
df_TransportationCosts = CSV.File(joinpath(@__DIR__, "TransportationCosts.csv")) |> DataFrame
df_TransportationDistance = CSV.File(joinpath(@__DIR__, "TransportationDistance.csv")) |> DataFrame
df_PricesPower = CSV.File(joinpath(@__DIR__, "PricesPower.csv")) |> DataFrame
df_PricesGas = CSV.File(joinpath(@__DIR__, "PricesGas.csv")) |> DataFrame
df_PricesBiomass = CSV.File(joinpath(@__DIR__, "PricesBiomass.csv")) |> DataFrame
df_PricesMethanol = CSV.File(joinpath(@__DIR__, "PricesMethanol.csv")) |> DataFrame
df_RegionsNames = CSV.File(joinpath(@__DIR__, "RegionsNames.csv")) |> DataFrame

df_name = [df_AnnualPowerDemand, df_AnnualHydrogenDemand, df_PowerDemandProfile, df_WindOnshoreTimeseries, df_SolarTimeseries, df_WindOffshoreTimeseries]
#replacing "," with "." in the dataframe and transforming it into Floats
for col in 3:size(df_AnnualPowerDemand)[2]
    for row in 1:size(df_AnnualPowerDemand)[1]
        df_AnnualPowerDemand[row, col] = replace.(df_AnnualPowerDemand[row, col], "," => ".")
    end
    transform!(df_AnnualPowerDemand, names(df_AnnualPowerDemand)[col] .=> ByRow(x -> parse(Float64, x)) .=> names(df_AnnualPowerDemand)[col])
end
for df in df_name[3:end]
    for col in 2:size(df)[2]
        if df[2, col] isa Number # if it is of type number, (some csv files store numbers already)
            continue
        else
            for row in 1:size(df)[1]
                df[row, col] = replace.(df[row, col], "," => ".")
            end
            transform!(df, names(df)[col] .=> ByRow(x -> parse(Float64, x)) .=> names(df)[col])
        end
    end
end
for col in 2:size(df_AnnualHydrogenDemand)[2]
    for row in 1:size(df_AnnualHydrogenDemand)[1]
        df_AnnualHydrogenDemand[row, col] = df_AnnualHydrogenDemand[row, col] * 1000000000 # convert PJ to MJ
        println(df_AnnualHydrogenDemand[row, col])
    end
end
timesteps = 1:8760
technologies = ["Electrolysis", "Gasification", "SMR", "POM"]
fuels = ["Power", "Gas", "Biomass", "Methanol"]
reg = 1:16

prod_power = Dict()

### NEU: PowerProd in MJ nach Bundesländern ###
STS = zeros(timesteps,reg)
WON = zeros(timesteps,reg)
WOF = zeros(timesteps,reg)
PowerProd = zeros(timesteps,reg)

for t in timesteps
    for r in reg
        STS[t,r] = df_CapacityInstalledGW[1,1+r]*df_SolarTimeseries[t,1+r]*3600000*2
        WON[t,r] = df_CapacityInstalledGW[2,1+r]*df_WindOnshoreTimeseries[t,1+r]*3600000*2
        WOF[t,r] = df_CapacityInstalledGW[3,1+r]*df_WindOffshoreTimeseries[t,1+r]*3600000*2
    end
end

for r in reg
    for i in timesteps
        PowerProd[i,r] = STS[i,r] + WON[i,r] + WOF[i,r]
    end
end

### NEU: Stromnachfrage nach Bundesländern ###

PowerDem = zeros(timesteps,16)
for r in reg
    for i in timesteps
        PowerDem[i,r] = df_AnnualPowerDemand[r,8] * df_PowerDemandProfile[i,r+1] * 1000000000
    end
end

### Set prices as x_fuels ###
x_fuels = Dict()
for i in timesteps
    for r in reg
        x_fuels["Power", i, r] = df_PricesPower[i, r+1]
        x_fuels["Gas", i, r] = df_PricesGas[i, r+1]
        x_fuels["Biomass", i, r] = df_PricesBiomass[i, r+1]
        x_fuels["Methanol", i, r] = df_PricesMethanol[i, r+1]
    end
end

#costs for technology e.g. personell costs, rent 
c_prod = Dict()
c_prod["Electrolysis"] = 0 
c_prod["Gasification"] = 0 
c_prod["SMR"] = 0.0016 # includes Carbon Tax
c_prod["POM"] = 0.00275

# Investment costs
c_invCostsCap = Dict()
c_invCostsCap["Electrolysis"] = 150000000 / 20 /240000 # euros to usd, 20 years depreciation
c_invCostsCap["Gasification"] = 3000000 / 20 / 5000
c_invCostsCap["SMR"] = 241000000 / 20 / 300000 # without Carbon capture
c_invCostsCap["POM"] = 100000000 / 20 / 100000


# how many MJ do we need for one MJ of hydrogen (efficiency)
input = Dict()
for j in fuels, i in technologies
    input[j, i] = 0
end
input["Power", "Electrolysis"] = 1.5
input["Power", "SMR"] = 0.0041
input["Gas", "SMR"] = 1.266
input["Biomass", "Gasification"] = 3.6
input["Methanol", "POM"] = 1.33

# capacity of a resource at any timestep (basically unlimited due to lacking reasonable alternatives)
r_max = Dict()
for j in fuels, t in timesteps
    r_max[j, t] = 10^17
end


# how many MegaJoule can be produced maximum in a single location
prod_max = Dict(zip(technologies, [240000, 5000, 300000, 100000]))

#Hydrogen Demand for regions and timesteps
hyd_reg_2040 = zeros((timesteps, 16))

for r in 1:16
    for t in timesteps
        hyd_reg_2040[t,r] = df_AnnualHydrogenDemandRegions[r+1,6]*df_PowerDemandProfile[t,r+1]*1000000000#*20 #(20*20TWh = 400 TWh)
    end
end

# costs for usage of external power 
c_external_power_prod = zeros(length(timesteps))
for timestep in timesteps
    #c_external_power_prod[timestep] = 0.9 - x_fuels["Power",timestep,1] # difference to current price for power, unit: USD/MJ
    c_external_power_prod[timestep] = x_fuels["Power",timestep,1]*10
end

#maximum storage per tank
st_max = 100000

#investmentcosts for one storage tank
c_inv_storage = (st_max * 70) / 20

#maximum storage per tank in MJ
st_max = 100000

#investmentcosts for one storage tank in USD
c_inv_storage = (st_max * 70) / 20 # 20 years depreciation and 70 per MJ storage
c_inv_storage = 70/20

#various efficiencies and losses
reg = 1:16
timesteps = 1:8760
storeEff = 1
storingEff = 1
transLoss = 0.00000263
transLossPow = 0.00008125
reElEff = 0.6
reElCosts = 0
compression = 0.12

#Transportation costs and transportation losses
TC = zeros(1:16, 1:16)
TD = zeros(1:16, 1:16)

for r in reg
    for i in reg
        TC[r,i] = df_TransportationCosts[r,i+1]
        TD[r,i] = df_TransportationDistance[r,i+1]
    end
end



#------- Model starts --------

ESM_H2_Base = Model(Cbc.Optimizer)


@variable(ESM_H2_Base, use_external_power[timesteps,reg] >= 0) # how many kwh of external power should be used
@variable(ESM_H2_Base, max_power[timesteps,reg] >= 0) # how much power can be used
@variable(ESM_H2_Base, prod[technologies, timesteps, reg] >= 0) #total hydrogen production in every hour
@variable(ESM_H2_Base, st[timesteps,reg] >= 0) #total storage level 
@variable(ESM_H2_Base, prodCapacity[technologies,reg] >= 0) #number of plants
@variable(ESM_H2_Base, capStorage[reg] >= 0) #number of storage tanks
@variable(ESM_H2_Base, rv[timesteps, reg] >= 0)
@variable(ESM_H2_Base, tr[timesteps,reg,reg] >= 0)
@variable(ESM_H2_Base, trp[timesteps,reg,reg] >= 0)


@objective(ESM_H2_Base, Min,
    sum(
        ((x_fuels[fuel,timestep,r] * input[fuel,technology] + c_prod[technology]) + x_fuels["Power",timestep,r]*compression) * prod[technology,timestep,r]
        for r in reg, fuel in fuels, technology in technologies, timestep in timesteps) + 
    sum(c_inv_storage * capStorage[r] for r in reg) + 
    sum(prodCapacity[technology,r] * c_invCostsCap[technology] for technology in technologies, r in reg) +
    sum(c_external_power_prod[timestep] * use_external_power[timestep,r] for timestep in timesteps, r in reg) + 
    sum(TC[r,re] * tr[t,r,re] for t in timesteps, r in reg, re in reg) #+
    #sum(rv[timestep,r] * reElCosts for timestep in timesteps, r in reg) #since reElCosts = 0, the objective function was shortened to save computational time
    )



#limit the production
@constraint(ESM_H2_Base, ProductionLimit[t in 1:length(timesteps), i in technologies, r in reg],
    prod[i, t, r] <= prodCapacity[i,r]
)

#limit the resources
@constraint(ESM_H2_Base, ResourceLimit[t in timesteps, j in fuels; !occursin(j, "Power")],
    sum(input[j, i] * prod[i, t, r] for i in technologies, r in reg) <= r_max[j, t]
)

#How much power can be used
@constraint(ESM_H2_Base, PowerLimit[hour in 2:length(timesteps), r in reg],
PowerProd[hour, r] - PowerDem[hour, r] + use_external_power[hour, r] + rv[hour, r]*reElEff + sum(trp[hour, re, r] * (1-transLossPow*TD[r,re]) for re in reg) >=
sum(input["Power", technology] * prod[technology, hour, r] for technology in technologies) + sum(prod[technology, hour, r]*compression for technology in technologies)
+ sum(trp[hour, r, re] for re in reg)
 ## power is in pj nut we need mj
)

#power balance, make max_power dependent of use_external_power
#@constraint(ESM_H2, PowerBilance[hour in 2:length(timesteps), r in reg],
#    max_power[hour, r] <= PowerProd[hour, r] - PowerDem[hour, r] + use_external_power[hour, r] + rv[hour, r]
#)


#Storage balance
@constraint(ESM_H2_Base, StorageBilance[t in 2:length(timesteps), r in reg],
    st[t, r] == st[t-1, r]*storeEff + sum(prod[i, t, r] for i in technologies) - hyd_reg_2040[t, r] - sum(tr[t,r,re] for re in reg) + sum(tr[t, re, r]*(1-TD[re,r]*transLoss) for re in reg) - rv[t, r]#*storingEff
)

#storage at beginning
@constraint(ESM_H2_Base, storageFirstPeriod[r in reg], 
    st[1,r] == 0.1*capStorage[r] #sum(prod[i, 1, r] for i in technologies) - hyd_reg_2035[1, r] - sum(tr[1,r,re] for re in reg) + sum(tr[1, re, r] for re in reg)
)

#storage am Ende
@constraint(ESM_H2_Base, storageLastPeriod[r in reg],
    st[length(timesteps),r] == st[1,r]
)

# storage capacity
@constraint(ESM_H2_Base, StorageLimit[t in 1:length(timesteps), r in reg],
    st[t,r] <= capStorage[r]
)

#@constraint(ESM_H2, StorageGgn[t in 1:length(timesteps), r in reg],
#    st[t,r] >= 0
#)

#@constraint(ESM_H2, Rückverstromung[t in 2:length(timesteps), r in reg],
#    rv[t,r] <= st[t-1,r] - hyd_reg_2035[t,r] + sum(prod[i, t, r] for i in technologies)
#)

optimize!(ESM_H2_Base)
objective_value(ESM_H2_Base)
objective_value(ESM_H2)
#-----------------preparing results for plotting and further analysis ------------------------
obj_neu = value.(objective_value(ESM_H2_Base))
obj_alt = value.(objective_value(ESM_H2))

#convert the dictionary into an array
prod_Electrolysis = zeros(length(timesteps), 1:16)
prod_Gasification = zeros(length(timesteps),1:16)
prod_SMR = zeros(length(timesteps),1:16)
prod_POM = zeros(length(timesteps),1:16)
storage = zeros(length(timesteps),reg)
r_treg = zeros(length(timesteps), reg)
external_power = zeros(length(timesteps),reg)
transportation = zeros(length(timesteps), reg, reg)
transpP = zeros(length(timesteps), reg, reg)

for t in timesteps
        println(t)
    for r in reg
        prod_Electrolysis[t,r] = value.(prod["Electrolysis",t,r])
        prod_Gasification[t,r] = value.(prod["Gasification",t,r])
        prod_SMR[t,r] = value.(prod["SMR",t,r])
        prod_POM[t,r] = value.(prod["POM",t,r])
        storage[t,r]=value.(st[t,r])
        r_treg[t,r]=value.(rv[t,r])
        external_power[t,r]=value.(use_external_power[t,r])
        for re in reg
            transportation[t,r,re]=value.(tr[t,r,re])
            #transpP[t,r,re] = value.(trp[t,r,re])
        end
    end
end

prod_total = hcat(prod_Electrolysis,prod_Gasification,prod_SMR,prod_POM)

prod_total_bl = zeros(reg,1:4)
for r in reg
    prod_total_bl[r,1] = sum(prod_Electrolysis[i,r] for i in timesteps)
    prod_total_bl[r,2] = sum(prod_Gasification[i,r] for i in timesteps)
    prod_total_bl[r,3] = sum(prod_SMR[i,r] for i in timesteps)
    prod_total_bl[r,4] = sum(prod_POM[i,r] for i in timesteps)
end