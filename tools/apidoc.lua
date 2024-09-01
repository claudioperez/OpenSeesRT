
-- pandoc's List type
local List = require 'pandoc.List'


--- Filter function for code blocks
function Div (cb)
  local apidoc = cb.attributes['apidoc']
  if not apidoc then
    return
  end
  io.stderr:write("apidoc.lua: " .. apidoc .. "\n")

  local blocks = List:new()
  fh = io.popen("python -m opensees.emit.apidoc " .. apidoc)
  local contents = pandoc.read(fh:read '*all', "markdown").blocks
  fh:close()
  -- contents = pandoc.walk_block(
  --   pandoc.Div(contents), { }
  -- ).content
  blocks:extend(contents)
  -- blocks:extend(cb)
  -- return blocks
  return blocks
end




  
