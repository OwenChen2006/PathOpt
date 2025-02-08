import javax.swing.*;
import java.awt.*;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.util.ArrayList;

public class CitySelector extends JPanel {
    private JTextField searchField;
    private JList<String> citiesList;
    private DefaultListModel<String> listModel;
    private ArrayList<String> allCities;
    private String selectedCity;
    private JDialog dialog;

    public CitySelector(ArrayList<String> cities) {
        this.allCities = cities;
        this.setLayout(new BorderLayout());
        
        // Initialize components
        searchField = new JTextField(20);
        listModel = new DefaultListModel<>();
        citiesList = new JList<>(listModel);
        citiesList.setVisibleRowCount(10);
        
        // Add all cities initially
        updateList("");
        
        // Add search functionality
        searchField.getDocument().addDocumentListener(new javax.swing.event.DocumentListener() {
            public void insertUpdate(javax.swing.event.DocumentEvent e) { updateList(searchField.getText()); }
            public void removeUpdate(javax.swing.event.DocumentEvent e) { updateList(searchField.getText()); }
            public void changedUpdate(javax.swing.event.DocumentEvent e) { updateList(searchField.getText()); }
        });
        
        // Add double-click functionality
        citiesList.addMouseListener(new MouseAdapter() {
            public void mouseClicked(MouseEvent e) {
                if (e.getClickCount() == 2 && citiesList.getSelectedValue() != null) {
                    selectedCity = citiesList.getSelectedValue();
                    if (dialog != null) {
                        dialog.dispose();
                    }
                }
            }
        });
        
        // Layout
        this.add(new JLabel("Search:"), BorderLayout.NORTH);
        this.add(searchField, BorderLayout.CENTER);
        this.add(new JScrollPane(citiesList), BorderLayout.SOUTH);
    }
    
    private void updateList(String searchText) {
        listModel.clear();
        String search = searchText.toLowerCase();
        for (String city : allCities) {
            if (city.toLowerCase().contains(search)) {
                listModel.addElement(city);
                if (listModel.size() >= 100) break; // Limit displayed results for performance
            }
        }
    }
    
    public String getSelectedCity() {
        return selectedCity;
    }
    
    // Static method to show dialog
    public static String showDialog(ArrayList<String> cities, String title) {
        JFrame parent = null;
        CitySelector selector = new CitySelector(cities);
        JDialog dialog = new JDialog(parent, title, true);
        selector.dialog = dialog;  // Store reference to dialog
        
        dialog.setContentPane(selector  );
        dialog.pack();
        dialog.setLocationRelativeTo(parent);
        dialog.setVisible(true);
        
        return selector.getSelectedCity();
    }
    
    // Example usage
    public static void main(String[] args) {
        // Sample data
        ArrayList<String> cities = new ArrayList<>();
        cities.add("New York, NY");
        cities.add("Los Angeles, CA");
        cities.add("Chicago, IL");
        // ... add your 30000 cities here
        
        SwingUtilities.invokeLater(() -> {
            String selected = CitySelector.showDialog(cities, null);
            if (selected != null) {
                System.out.println("Selected: " + selected);
            }
        });
    }
}