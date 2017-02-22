function geneticTSP()
%    2017-02-21 by  Tungpo  in 834 
    cityNum = 15;
    cityList = rand(2,cityNum);

    %参数
    maxGeneration = 2000;
    population = 100;
    pcross = 0.2;
    pmutant = 0.5;

    %初始种群
    currentGeneration = zeros(population,cityNum);
    for i = 1:population
        currentGeneration(i,:) = randperm(cityNum);
    end

    historyBest = zeros(maxGeneration,cityNum+1);
    nextGeneration = currentGeneration;
    for gn = 1:maxGeneration
        currentGeneration = nextGeneration;
    	%计算是适应度
        fitness = fitnessCal(currentGeneration);
        %存储历代最优个体
        [M,I] = max(fitness);
        historyBest(gn,:) = [currentGeneration(I,:) M];
        %选择下一代个体
        nextGeneration = select(currentGeneration, fitness);
        %染色体交叉
        nextGeneration = crossover(nextGeneration);
        %染色体变异
        nextGeneration = mutant(nextGeneration);
        nextGeneration(1,:) = historyBest(gn,1:cityNum);
        %显示目前为止最优解
        drawPath(historyBest,gn);
    end

    function fitness = fitnessCal(currentG)
        [popu, ~] = size(currentG);
        fitness = zeros(popu,1);
        for indiv = 1:popu
            fitness(indiv,1) = cal(currentG(indiv,:));
        end
    end



    function fit = cal(que)
        [~,cnum] = size(que);
        dis = 0;
        for ii = 2:cnum
            dis = dis + sqrt( (cityList(1,que(1,ii-1)) - cityList(1,que(1,ii)))* (cityList(1,que(1,ii-1)) - cityList(1,que(1,ii)))...
                                +(cityList(2,que(1,ii-1)) - cityList(2,que(1,ii)))*(cityList(2,que(1,ii-1)) - cityList(2,que(1,ii))));
        end
        dis = dis + sqrt( (cityList(1,que(1,1)) - cityList(1,que(1,cityNum)))* (cityList(1,que(1,1)) - cityList(1,que(1,cityNum)))...
                                +(cityList(2,que(1,1)) - cityList(2,que(1,cityNum)))*(cityList(2,que(1,1)) - cityList(2,que(1,cityNum))));
        fit = 1000/dis;
    end
    
    function drawPath(hist,gn)
        [Maxfit, maxnum]=max(hist(1:gn,cityNum+1));
        que = hist(maxnum,:);
        disp(['当前代数：' num2str(gn) '当前适应度：' num2str(Maxfit)]);
        for ii = 2:cityNum
            plot([cityList(1,que(1,ii-1)) cityList(1,que(1,ii))],...
                 [cityList(2,que(1,ii-1)) cityList(2,que(1,ii))],...
                 'Marker','.',...
                 'color','r');
             hold on;
        end
        plot(   [cityList(1,que(1,1)) cityList(1,que(1,cityNum))],...
                 [cityList(2,que(1,1)) cityList(2,que(1,cityNum))],...
                 'Marker','.',...
                 'color','r');
        hold off;
        pause(0.05);     
        
    end

    
    function newGeneration = select(oldGeneration, fits)
        [popu,cn] =size(oldGeneration);
        newGeneration = zeros(popu,cn);
        %计算概率
        pselect = fits(:,1)/sum(fits(:,1));
        psum = zeros(popu+1 , 1);
        for ii = 2:popu+1
            psum(ii,1) = psum(ii-1,1) + pselect(ii-1,1);
        end
        %选择后代
        for ii = 1:popu
            seed = rand(1);
            for jj = 1:popu+1
                if seed < psum(jj,1)
                    newGeneration(ii,:) = oldGeneration(jj-1,:);
                    break;
                end
            end
        end
        
    end
    
    function newGeneration = crossover(oldGeneration)
        [popu, cnum] = size(oldGeneration);
        mark = floor(popu/2);
        newGeneration = oldGeneration;
        p = 0;
        %交叉
        for ii = 1:mark
            p = rand(1);
            if p < pcross
                location = round(rand(1,2)*cnum);
                sort(location);
                for jj = 1:cnum
                    if jj < location(1,1)|| jj > location(1,2)
                        newGeneration(ii,jj) = oldGeneration(ii, jj);
                        newGeneration(ii+mark,jj) = oldGeneration(ii+mark, jj);
                    else
                        newGeneration(ii, jj) = oldGeneration(ii+mark, jj);
                        newGeneration(ii+mark, jj) = oldGeneration(ii, jj);
                    end
                end
            end
        end
        %检查合法性
        for ii = 1:popu
            chrom = newGeneration(ii,:);
            num = zeros(1,cnum);
            buffer = zeros(1,cnum);
            buffernum = 1;
            %统计出现次数
            for jj = 1:cnum
                num(1,chrom(1,jj)) = num(1,chrom(1,jj)) + 1;
            end
            %存储未出现过的数字
            for jj = 1:cnum
               if num(1,jj) == 0
                   buffer(1,buffernum) = jj;
                   buffernum = buffernum + 1;
               end
            end
            %找到未出现过数字目前可以出现的位置
            for jj = 1:cnum
                if num(1,chrom(1,jj)) > 1
                    num(1,chrom(1,jj)) = num(1,chrom(1,jj)) - 1;
                    chrom(1,jj) = 0;
                end
            end
            %调整:将未出现过的数字填入目前为0的位置
            buffernum = 1;
            for jj = 1:cnum
                if chrom(1,jj) == 0
                    chrom(1,jj) = buffer(1,buffernum);
                    buffernum = buffernum + 1;
                end
            end
            newGeneration(ii,:) = chrom;
        end
        
    end


    function newGeneration = mutant(oldGeneration)
        newGeneration = oldGeneration;
        [popu, cnum] = size(oldGeneration);
        d = 0;
        for ii = 1:popu
            d=rand(1);
            if d < pmutant
                time = round(rand(1)*cnum);
                for jj = 1:time
                    location = ceil(rand(1,2)*cnum);
                    temp = newGeneration(ii,location(1,1));
                    newGeneration(ii,location(1,1)) = newGeneration(ii,location(1,2));
                    newGeneration(ii,location(1,2)) = temp;
                end
            end
        end
    end

end







