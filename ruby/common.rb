#
module Common
  # Load the common.sh data in a hash
  def load_common
    h={}
    common_file = File.dirname(__FILE__) + "/../bash/common.sh"
    File.open(common_file).each_line do |l| 
      if l =~ /^[\w_-]+=["-\/\$\w_.]+\"$/
        var, val = l.split('=')
        h[var] = val.gsub(/"/, '').chomp
        #puts var + ":" + h[var]
      end
    end 
    h
  end

  # check for data in STDIN
  def data_in_stdin?
    require 'fcntl'
    STDIN.fcntl(Fcntl::F_GETFL, 0) == 0 ? true : false
  end

  def log(msg, same_line=0, c_return=1)
    cr = c_return  == 1 ? "\n" : ""
    sl = same_line == 1 ? "\r" : ""
    $stderr.printf "#{sl}%s#{cr}", msg
  end
end

module Help
  @usage_text = ""

  def self.set_usage_text(ut)
    @usage_text = ut
  end

  def self.error(msg)
    puts "ERROR: #{msg}"
    usage
    exit 1
  end

  def self.usage
    puts @usage_text
  end
end

module Statistics
  def variance(population) 
    n = 0         
    mean = 0.0  
    s = 0.0  
    population.each do |x|  
      n = n + 1     
      delta = x - mean
      mean = mean + (delta / n)    
      s = s + delta * (x - mean)  
    end

    # if you want to calculate std deviation
    # of a sample change this to "s / (n-1)"
    return s / n                          
  end

  # calculate the standard deviation of a population
  # accepts: an array, the population       
  # returns: the standard deviation    
  def standard_deviation(population)  
    Math.sqrt(variance(population))  
  end             

  def mean(a)     
    sum  = a.inject(0) {|r,i| r.to_i + i.to_i }
    sum.to_f / a.size.to_f  
  end 
end

