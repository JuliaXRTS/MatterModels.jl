"""
    GridLookupTable(
        cols::T1,
        rows::T2,
        values::T3
    ) where {T1<:AbstractVector{T},T2<:AbstractVector{T},T3<:AbstractMatrix{T}} where {T<:Real}

TBW
"""
struct GridLookupTable{T1, T2, T3}
    rows::T1
    cols::T2
    values::T3

    function GridLookupTable(
            cols::T1, rows::T2, values::T3
        ) where {
            T1 <: AbstractVector{T}, T2 <: AbstractVector{T}, T3 <: AbstractMatrix{T},
        } where {T <: Real}
        len_arg_rows = length(rows)
        len_arg_cols = length(cols)
        len_val_rows, len_val_cols = size(values)

        len_arg_cols == len_val_cols || throw(
            ArgumentError(
                "Length of column arguments <$len_arg_cols> does not match the the first dimension of the values matrix <$len_val_cols>.",
            ),
        )

        len_arg_rows == len_val_rows || throw(
            ArgumentError(
                "Length of row arguments <$len_arg_rows> does not match the the second dimension of the values matrix <$len_val_rows>.",
            ),
        )

        # TODO: implement check, if rows and cols are in ascending order

        return new{T1, T2, T3}(rows, cols, values)
    end
end

# TODO:
# - change signature of lookup to lookup(GLT, method, column,row)
# - check behaviour for sampling outside of the lookup
# - add some I/O for lookup tables (simple saving and loading)

########
# extrapolation methods
########
abstract type AbstractLookupMethod end
struct InterpolExtrapol <: AbstractLookupMethod end
struct InterpolEndValue <: AbstractLookupMethod end
struct NearestInput <: AbstractLookupMethod end
struct BelowInput <: AbstractLookupMethod end
struct AboveInput <: AbstractLookupMethod end

"""
    lookup(
    method::InterpolExtrapol,
    column::Real,
    row::Real,
    table::GridLookupTable)

TBW
"""
function lookup(method::InterpolExtrapol, column::Real, row::Real, table::GridLookupTable)
    # Binary search in rows vector
    rowlow = searchsortedfirst(table.rows, row)
    columnlow = searchsortedfirst(table.cols, column)


    # Behaviour according to method
    if rowlow >= length(table.rows)
        rowlow = length(table.rows) - 1
    end
    if columnlow >= length(table.cols)
        columnlow = length(table.cols) - 1
    end

    rowhigh = rowlow + 1
    columnhigh = columnlow + 1
    # Calculate column and row factors
    rowfactor = (row - table.rows[rowlow]) / (table.rows[rowhigh] - table.rows[rowlow])
    columnfactor =
        (column - table.cols[columnlow]) / (table.cols[columnhigh] - table.cols[columnlow])
    # For row low
    outputcolumnΔlow = table.values[rowlow, columnhigh] - table.values[rowlow, columnlow]
    # For row high
    outputcolumnΔhigh = table.values[rowhigh, columnhigh] - table.values[rowhigh, columnlow]
    # For this column
    outputrowΔ =
        table.values[rowhigh, columnlow] + outputcolumnΔhigh * columnfactor -
        (table.values[rowlow, columnlow] + outputcolumnΔlow * columnfactor)
    return table.values[rowlow, columnlow] +
        columnfactor * outputcolumnΔlow +
        rowfactor * outputrowΔ
end

"""
    lookup(
    method::InterpolEndValue,
    column::Real,
    row::Real,
    table::GridLookupTable)

TBW
"""
function lookup(method::InterpolEndValue, column::Real, row::Real, table::GridLookupTable)
    # Binary search in rows vector
    rowlow = searchsortedfirst(table.rows, row)
    columnlow = searchsortedfirst(table.cols, column)

    if rowlow == length(table.rows)
        rowlow -= 1
    end
    if columnlow == length(table.cols)
        columnlow -= 1
    end
    rowhigh = rowlow + 1
    columnhigh = columnlow + 1
    if column >= table.cols[columnhigh]
        columnfactor = 0
        outputcolumnΔlow = 0
        columnlow = columnhigh
    elseif column <= table.cols[columnlow]
        columnfactor = 0
        outputcolumnΔlow = 0
    else
        # For row low
        outputcolumnΔlow =
            table.values[rowlow, columnhigh] - table.values[rowlow, columnlow]
    end
    if row >= table.rows[rowhigh]
        rowfactor = 0
        outputrowΔ = 0
        rowlow = rowhigh
    elseif row <= table.rows[rowlow]
        rowfactor = 0
        outputrowΔ = 0
    else
        # For row high
        outputcolumnΔhigh =
            table.values[rowhigh, columnhigh] - table.values[rowhigh, columnlow]
        # Calculate row factor
        rowfactor = (row - table.rows[rowlow]) / (table.rows[rowhigh] - table.rows[rowlow])
        # For this column
        outputrowΔ =
            table.values[rowhigh, columnlow] + outputcolumnΔhigh * columnfactor -
            (table.values[rowlow, columnlow] + outputcolumnΔlow * columnfactor)
    end
    return table.values[rowlow, columnlow] +
        columnfactor * outputcolumnΔlow +
        rowfactor * outputrowΔ
end

"""
    lookup(
    method::NearestInput,
    column::Real,
    row::Real,
    table::GridLookupTable)

TBW
"""
function lookup(method::NearestInput, column::Real, row::Real, table::GridLookupTable)
    # Binary search in rows vector
    rowlow = searchsortedfirst(table.rows, row)
    columnlow = searchsortedfirst(table.cols, column)
    if rowlow == length(table.rows)
        rowlow -= 1
    end
    if columnlow == length(table.cols)
        columnlow -= 1
    end
    rowhigh = rowlow + 1
    columnhigh = columnlow + 1
    if row < (table.rows[rowlow] + table.rows[rowhigh]) / 2
        nearestrow = rowlow
    else
        nearestrow = rowhigh
    end
    if column < (table.cols[columnlow] + table.cols[columnhigh]) / 2
        nearestcolumn = columnlow
    else
        nearestcolumn = columnhigh
    end
    return table.values[nearestrow, nearestcolumn]
end

"""
    lookup(
    method::BelowInput,
    column::Real,
    row::Real,
    table::GridLookupTable)

TBW
"""
function lookup(method::BelowInput, column::Real, row::Real, table::GridLookupTable)
    # Binary search in rows vector
    rowlow = searchsortedfirst(table.rows, row)
    columnlow = searchsortedfirst(table.cols, column)
    return table.values[rowlow, columnlow]
end

"""
    lookup(
    method::AboveInput,
    column::Real,
    row::Real,
    table::GridLookupTable)

TBW
"""
function lookup(method::AboveInput, column::Real, row::Real, table::GridLookupTable)
    # Binary search in rows vector
    rowlow = searchsortedfirst(table.rows, row)
    columnlow = searchsortedfirst(table.cols, column)
    if table.rows[rowlow] == row
        rowhigh = rowlow
    else
        rowhigh = rowlow + 1
    end
    if table.cols[columnlow] == column
        columnhigh = columnlow
    else
        columnhigh = columnlow + 1
    end
    return table.values[rowhigh, columnhigh]
end
