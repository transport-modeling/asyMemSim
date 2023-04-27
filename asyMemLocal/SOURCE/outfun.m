function stop = outfun(x,optimValues,state)
stop = false;
 
   switch state
       case 'init'
           hold on
       case 'iter'

           plot(optimValues.iteration,norm(optimValues.fval)^2,'o');
           text(optimValues.iteration+0.015,norm(optimValues.fval)^2+0.005,num2str(optimValues.iteration));

       case 'done'
           hold off
       otherwise
   end
end

