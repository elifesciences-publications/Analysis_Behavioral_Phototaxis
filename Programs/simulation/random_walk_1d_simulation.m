function[x2_max] =  random_walk_1d_simulation ( step_num, walk_num )

%*****************************************************************************80
%
%% RANDOM_WALK_1D_SIMULATION simulates a random walk in 1D.
%
%  Discussion:
%
%    The expectation should be that, the average distance squared X^2 
%    is equal to the time, or number of steps N.
%
%    Or, equivalently
%
%      average ( |X| ) = sqrt ( N )
%
%    The program makes a plot of both the average and the maximum values
%    of X^2 versus time.  The maximum value grows much more quickly,
%    and that curve is much more jagged.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    04 November 2009
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer STEP_NUM, the number of steps to take in one test.
%
%    Input, integer WALK_NUM, the number of times to repeat the walk.
%
  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'RANDOM_WALK_1D_SIMULATION:\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Simulate WALK_NUM random walks of STEP_NUM steps each.\n' );

  if ( nargin < 1 )
    step_num = input ( 'Enter STEP_NUM, the number of steps to take:  ' );
  elseif ( ischar ( step_num ) )
    step_num = str2num ( step_num );
  end

  if ( nargin < 2 )
    walk_num = input ( 'Enter WALK_NUM, the number of walks to make:  ' );
  elseif ( ischar ( walk_num ) )
    walk_num = str2num ( walk_num );
  end

  time = 0 : step_num;
  x2_ave = zeros(step_num+1,1);
  x2_max = zeros(step_num+1,1);
%
%  Take WALK_NUM walks, each of STEP_NUM random +1 or -1 steps.
%
  for walk = 1 : walk_num

    x = 0;

    for step = 1 : step_num

      r = rand ( );
%
%  Take the stap.
%
      if ( r <= 0.5 )
        x = x - 1;
      else
        x = x + 1;
      end
%
%  Update the average and max.
%
      x2_ave(step+1) = x2_ave(step+1) + x^2;
      x2_max(step+1) = max ( x2_max(step+1), x^2 );

    end

  end
%
%  Average over the number of walks.
%
  x2_ave(:,:) = x2_ave(:,:) / walk_num;
%
%  Plot the results.
%
  plot ( time, x2_ave, time, x2_max, 'LineWidth', 2 );

  title ( '1D Random Walk - Max and average of distance^2 versus time' )
  xlabel ( 'Time' )
  ylabel ( 'Distance Squared' )
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'RANDOM_WALK_1D_SIMULATION:\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  return
end
function timestamp ( )

%*****************************************************************************80
%
%% TIMESTAMP prints the current YMDHMS date as a timestamp.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    14 February 2003
%
%  Author:
%
%    John Burkardt
%
  t = now;
  c = datevec ( t );
  s = datestr ( c, 0 );
  fprintf ( 1, '%s\n', s );

  return
end