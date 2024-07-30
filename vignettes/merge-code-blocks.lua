local code_block_sep = '\n'
function Blocks (blocks)
for i = #blocks, 2, -1 do -- start at end of list
  -- Both blocks must be code blocks and share the same primary class
if blocks[i - 1].t == 'CodeBlock' and
blocks[i].t == 'CodeBlock' and
blocks[i - 1].classes[1] == blocks[i].classes[1] then
blocks[i - 1].text = blocks[i - 1].text ..
code_block_sep ..
blocks[i].text
blocks:remove(i)
end
end
return blocks
end
