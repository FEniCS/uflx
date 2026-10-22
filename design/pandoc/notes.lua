-- Adapt the design notes for HTML and PDF output.
local ext = "html"
local note = ""

function Meta(meta)
  if meta["link-ext"] then
    ext = pandoc.utils.stringify(meta["link-ext"])
  end
  if meta["note"] then
    note = pandoc.utils.stringify(meta["note"])
  end
end

function Link(el)
  local target = el.target
  if target:match("^%a+://") then
    return nil
  end
  local path, anchor = target:match("^([^#]*%.md)(#?.*)$")
  if path then
    if ext == "pdf" and path == note and anchor ~= "" then
      -- A link inside one PDF must be a bare anchor to resolve.
      el.target = anchor
    else
      el.target = path:gsub("%.md$", "." .. ext) .. anchor
    end
  end
  return el
end

-- Wide tables need column widths, or LaTeX will not wrap the cells.
function Table(tbl)
  if ext ~= "pdf" then
    return nil
  end
  local n = #tbl.colspecs
  local longest = {}
  for i = 1, n do
    longest[i] = 18
  end
  local function measure(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if i <= n then
          local len = #pandoc.utils.stringify(cell.contents)
          longest[i] = math.max(longest[i], math.min(len, 80))
        end
      end
    end
  end
  measure(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do
    measure(body.body)
  end
  local total = 0
  for i = 1, n do
    total = total + longest[i]
  end
  if total < 90 then
    return nil
  end
  for i = 1, n do
    tbl.colspecs[i] = { tbl.colspecs[i][1], longest[i] / total }
  end
  return tbl
end

-- Long inline code cannot break in LaTeX, so offer breaks after separators.
function Code(el)
  if ext ~= "pdf" or #el.text < 24 then
    return nil
  end
  local out = pandoc.List()
  for piece in el.text:gmatch("[^,%(%[ ]*[,%(%[ ]?") do
    if piece ~= "" then
      out:insert(pandoc.Code(piece, el.attr))
      out:insert(pandoc.RawInline("latex", "\\allowbreak{}"))
    end
  end
  return out
end

-- Blocks without a language skip fancyvrb, which does the line breaking.
function CodeBlock(el)
  if ext == "pdf" and #el.classes == 0 then
    el.classes = { "default" }
    return el
  end
end

return { { Meta = Meta }, { Link = Link, Table = Table, Code = Code, CodeBlock = CodeBlock } }
