# TODO
# - rename the type (adopt MCSS)
# - build a better interface of data informed matter models
# - update lookup function (signature)
struct DataDrivenMatterModel <: AbstractMatterModel
    lt::GridLookupTable
    method::AbstractLookupMethod
end
lookup_table(dsf::DataDrivenMatterModel) = dsf.lt
lookup_method(dsf::DataDrivenMatterModel) = dsf.method

function DataDrivenMatterModel(datadict::Dict, method::AbstractLookupMethod = InterpolExtrapol())
    lt = GridLookupTable(
        datadict["omega_me"],
        datadict["q_me"],
        datadict["dsf"]',
    )
    return DataDrivenMatterModel(lt, method)
end

function dynamic_structure_factor(
        sys::DataDrivenMatterModel,
        om_q::NTuple{2, T},
    ) where {T <: Real}
    om, q = @inbounds om_q

    return lookup(lookup_method(sys), om, q, lookup_table(sys))
end
