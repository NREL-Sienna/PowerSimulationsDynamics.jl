get_inner_vars(dynamic_device::PSY.DynamicInjection) = dynamic_device.ext[INNER_VARS]

function _get_value(v::ForwardDiff.Dual)
    return v.value
end

function set_inner_vars(dynamic_device, VAR, value::Float64)
    get_inner_vars(dynamic_device)[VAR] = value
end

function set_inner_vars(dynamic_device, VAR, value)
    set_inner_vars(dynamic_device, VAR, _get_value(value))
end
