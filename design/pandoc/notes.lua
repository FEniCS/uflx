-- Adapt the design notes for PDF output.
local note = ""

function Meta(meta)
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
    if path == note and anchor ~= "" then
      -- A link inside one PDF must be a bare anchor to resolve.
      el.target = anchor
    else
      el.target = path:gsub("%.md$", ".pdf") .. anchor
    end
  end
  return el
end

-- Wide tables need column widths, or LaTeX will not wrap the cells.
-- Each column gets room for its longest word, and the rest by content.
local LINE = 100

function Table(tbl)
  local n = #tbl.colspecs
  local longest, word = {}, {}
  for i = 1, n do
    longest[i], word[i] = 4, 4
  end
  local function measure(rows)
    for _, row in ipairs(rows) do
      for i, cell in ipairs(row.cells) do
        if i <= n then
          local text = pandoc.utils.stringify(cell.contents)
          longest[i] = math.max(longest[i], math.min(#text, 80))
          for w in text:gmatch("%S+") do
            word[i] = math.max(word[i], math.min(#w + 2, 40))
          end
        end
      end
    end
  end
  measure(tbl.head.rows)
  for _, body in ipairs(tbl.bodies) do
    measure(body.body)
  end
  local total, base, extra = 0, 0, 0
  for i = 1, n do
    total = total + longest[i]
    base = base + word[i]
    extra = extra + (longest[i] - math.min(word[i], longest[i]))
  end
  if total < 90 then
    return nil
  end
  local spare = math.max(LINE - base, 0)
  local widths, sum = {}, 0
  for i = 1, n do
    local grow = extra > 0 and (longest[i] - math.min(word[i], longest[i])) / extra or 0
    widths[i] = word[i] + spare * grow
    sum = sum + widths[i]
  end
  for i = 1, n do
    tbl.colspecs[i] = { tbl.colspecs[i][1], widths[i] / sum }
  end
  return tbl
end

-- Long inline code cannot break in LaTeX, so offer breaks after separators.
function Code(el)
  if #el.text < 24 then
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
  if #el.classes == 0 then
    el.classes = { "default" }
    return el
  end
end

return { { Meta = Meta }, { Link = Link, Table = Table, Code = Code, CodeBlock = CodeBlock } }
