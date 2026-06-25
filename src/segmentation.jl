"""
    segment_views(profile, scale; overlap=false, bidirectional=true)

Split `profile` into views of length `scale`.

With `overlap=false` the profile is divided into adjacent non-overlapping
segments; when `bidirectional=true` a second pass starting from the end of the
profile is added so that the trailing samples left over by integer division are
also covered. With `overlap=true` a sliding window of unit step is used and
`bidirectional` is ignored.
"""
function segment_views(
        profile::AbstractVector, scale::Integer; overlap::Bool = false, bidirectional::Bool = true
    )
    profile_length = length(profile)
    scale >= 2 || throw(ArgumentError("scale must be at least 2"))
    scale <= profile_length ||
        throw(ArgumentError("scale must not exceed the profile length"))

    segments = typeof(view(profile, 1:scale))[]
    if overlap
        for start in 1:(profile_length - scale + 1)
            push!(segments, view(profile, start:(start + scale - 1)))
        end
    else
        segment_count = fld(profile_length, scale)
        for index in 0:(segment_count - 1)
            start = index * scale + 1
            push!(segments, view(profile, start:(start + scale - 1)))
        end
        if bidirectional
            offset = profile_length - segment_count * scale
            for index in 0:(segment_count - 1)
                start = offset + index * scale + 1
                push!(segments, view(profile, start:(start + scale - 1)))
            end
        end
    end

    isempty(segments) && throw(ArgumentError("no segments produced for scale $scale"))
    return segments
end
