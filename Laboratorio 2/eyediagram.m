function [] = eyediagram(outputpower, numb, x)

  n = ceil(numb/2);
  inter = floor(numel(outputpower)/n);

  i = 1;
  for j = 1 : inter : numel(outputpower)
    if (inter + j < numel(outputpower))
      pause
      figure(x)
      plot([1 : 1 : inter + 1], outputpower(j : inter + j))
      hold on
      i = i + 1
    end
  end
  figure(x)
  grid on

end

%FIXME partire dalla transizione forse, non da dove inizia
