#
module Common
  # Load the common.sh data in a hash
  def self.load_common
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
end
