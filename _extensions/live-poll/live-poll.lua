-- {{< live-poll M3.1 >}}
-- Embeds the live poll presenter view (poll/present.html) for one moment.
-- Optional: height="520px"
-- Questions for each moment are written in the live-poll Google Sheet.

-- "../" repeated once per folder between this document and the project root,
-- so the link works whatever URL the site is served from.
local function path_to_root()
  local input = quarto.doc.input_file
  local root = quarto.project.directory
  if not input or not root then return "" end
  local rel = input:sub(#root + 2)
  local _, depth = rel:gsub("/", "")
  return string.rep("../", depth)
end

return {
  ["live-poll"] = function(args, kwargs)
    local moment = pandoc.utils.stringify(args[1] or "")
    if moment == "" then
      return pandoc.Strong(pandoc.Str("live-poll: missing moment, e.g. {{< live-poll M3.1 >}}"))
    end
    local height = pandoc.utils.stringify(kwargs["height"] or "")
    if height == "" then height = "540px" end
    local src = path_to_root() .. "poll/present.html?m=" .. moment

    if quarto.doc.is_format("revealjs") then
      -- data-src: reveal.js loads the iframe when its slide is shown and
      -- unloads it when you move on, so hidden polls don't keep refreshing.
      return pandoc.RawBlock("html", string.format(
        '<iframe class="live-poll" data-src="%s" style="width:100%%;height:%s;border:0;" allow="fullscreen"></iframe>',
        src, height))
    elseif quarto.doc.is_format("html") then
      return pandoc.RawBlock("html", string.format(
        '<iframe class="live-poll" src="%s" style="width:100%%;height:%s;border:0;" allow="fullscreen" loading="lazy"></iframe>',
        src, height))
    end
    return pandoc.Para({ pandoc.Link("Open the live poll (" .. moment .. ")", src) })
  end
}
