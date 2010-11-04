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


